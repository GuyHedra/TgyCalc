import math
import olof_1 as o1
# from numpy import np.array, around, set_printoptions, matmul
import numpy as np
import matplotlib.pyplot as plt
import numbers
from scipy.spatial.transform import Rotation as R

# todo: consider using numpy arrays for all vectors and coordinates, hopefully can use their funcs instead of homegrown
""" TgyTools.py is a simplified version of polyTools intended for use within TgyCalc"""

stable_tol = 0.0001

def vector_add(v0, v1):
    return [c0 + c1 for c0, c1 in zip(v0, v1)]


def vector_mag(vector):
    """ return magnitude of a vector"""
    return sum([c ** 2 for c in vector]) ** 0.5


def vector_scalar_multiply(vector, scalar):
    # if isinstance(vector, list) and isinstance(scalar, numbers.Number):
        return [element * scalar for element in vector]
    # else:
    #     raise Exception('vector must be a list and scalar must be a number')


def dot_product(v0, v1):
    """returns the dot product of two vectors"""
    return sum([v0_coord * v1_coord for v0_coord, v1_coord in zip(v0, v1)])


def vec_cos(vec0, vec1):
    return math.acos(dot_product(vec0, vec1) / (vector_mag(vec0) * vector_mag(vec1)))


def cross_product(vector0, vector1):
    """Returns the cross product between two free vectors"""
    return [vector0[1] * vector1[2] - vector0[2] * vector1[1],
            vector0[2] * vector1[0] - vector0[0] * vector1[2],
            vector0[0] * vector1[1] - vector0[1] * vector1[0]]


def vec_angle(vec0, vec1):
    """ returns the acute angle between two free vectors """
    value = math.acos(dot_product(vec0, vec1) / (vector_mag(vec0) * vector_mag(vec1)))
    # if the angle is greater than pi / 2 then it is the obtuse angle, and we want the acute angle
    if value > math.pi / 2:
        value = math.pi - value
    return value

def normalize(vector):
    magnitude = sum([element ** 2 for element in vector]) ** 0.5
    if magnitude != 0:
        return [element / magnitude for element in vector]
    else:
        return vector


def vector_projection(projected_vector, unit_vector):
    """ return the projection of projected_vector onto unit_vector"""
    cos_alpha = dot_product(projected_vector, unit_vector) / (vector_mag(projected_vector) * vector_mag(unit_vector))
    return vector_scalar_multiply(normalize(unit_vector), cos_alpha * vector_mag(projected_vector))


def cartesian_coordinates(cyl_coords):
    """ cylindrical are [rho(r), phi(azimuth), z]"""
    rho, phi, z = cyl_coords
    return [rho * math.cos(phi), rho * math.sin(phi), z]


class Vertex:
    def __init__(self, cartesian=None, cylindrical=None, step=0.1):
        if cylindrical:
            self.coordinates = cartesian_coordinates(cylindrical)
        elif cartesian:
            self.coordinates = cartesian
        else:
            raise Exception('Either cartesian or cylindrical must be passed')
        self.strut = None
        self.tendons = []
        self.step = step

    def set_coordinates(self, coordinates, cyl_coords=None):
        if cyl_coords:
            self.coordinates = cartesian_coordinates(cyl_coords)
        else:
            self.coordinates = coordinates

    def set_strut(self, strut):
        self.strut = strut

    def add_tendon(self, tendon):
        self.tendons.append(tendon)

    def set_step(self, step):
        self.step = step

    def rotate(self, center, axis, angle):
        """ Rotates self about center and axis by angle radians.
        Translates self such that center is at the origin, rotates self by angle radians about axis and then
        translates self such that center is back at its original location"""
        point = self - center
        r = R.from_rotvec(angle * np.array(axis))
        self.coordinates = r.apply(point) + center.coordinates

    @property
    def cyl_coordinates_rad(self):
        """ returns [rho(r), phi(azimuth), z] """
        x = self.coordinates[0]
        y = self.coordinates[1]
        z = self.coordinates[2]
        r = (x ** 2 + y ** 2) ** 0.5
        theta = math.atan2(y, x)
        return [r, theta, z]

    @property
    def cyl_coordinates_deg(self):
        """ returns [rho(r), phi(azimuth), z] """
        x = self.coordinates[0]
        y = self.coordinates[1]
        z = self.coordinates[2]
        r = (x ** 2 + y ** 2) ** 0.5
        theta = math.degrees(math.atan2(y, x))
        return [r, theta, z]

    @property
    def f_is_lateral(self):
        return abs(dot_product(normalize(self.f_vector), normalize(self.strut.axis_vector(self)))) < 2 ** 0.5 / 2

    @property
    def f_is_longitudinal(self):
        return not self.f_is_lateral

    @property
    def members(self):
        return [self.strut] + self.tendons

    @property
    def f_vector(self):
        """ returns the vector sum af all tendon spring forces acting on this vertex """
        force_vector = [0, 0, 0]
        for tendon in self.tendons:
            force_vector = vector_add(tendon.f_vector(self), force_vector)
        return force_vector

    def distance(self, other):
        return sum(
            [(p0_coord - p1_coord) ** 2 for p0_coord, p1_coord in zip(self.coordinates, other.coordinates)]) ** 0.5

    def __sub__(self, other):
        if isinstance(other, Vertex):
            return [self_coord - other_coord for self_coord, other_coord in zip(self.coordinates, other.coordinates)]
        else:
            return [self_coord - other_coord for self_coord, other_coord in zip(self.coordinates, other)]


class Member:
    """ base class for tendons and strut"""

    def __init__(self, vertices):
        self.vertices = vertices
        # self.force = None
        self.xalglib_force = None

    def xalglib_force_vector(self, vertex):
        return [element * self.xalglib_force for element in self.unit_force_vector(vertex)]

    def set_xalglib_force(self, force):
        self.xalglib_force = force

    @property
    def current_length(self):
        return self.vertices[0].distance(self.vertices[1])

    @property
    def midpoint(self):
        """ Returns the midpoint of self"""
        return [(c0 + c1) / 2 for c0, c1 in zip(self.vertices[0].coordinates, self.vertices[1].coordinates)]

    def rotate(self, center_vertex, axis, angle):
        for vertex in self.vertices:
            vertex.rotate(center=center_vertex.coordinates, axis=axis, angle=angle)

    def axis_vector(self, initial_vertex):
        """ return a free vector parallel to the strut that terminates at the terminal_vertex arg"""
        # todo: move this func from strut to member
        if self.vertices[1] is initial_vertex:
            return normalize(self.vertices[0] - self.vertices[1])
        elif self.vertices[0] is initial_vertex:
            return normalize(self.vertices[1] - self.vertices[0])
        else:
            raise Exception('strut does not contain initial_vertex')


class Strut(Member):
    """ Compression member. Has no spring force and exerts no force."""

    def __init__(self, vertices):
        Member.__init__(self, vertices)
        self.member_type = 'Strut '

    def translate(self, step):
        """ sum vertex spring vectors and translate strut by step in said vector direction"""
        displacement = vector_scalar_multiply(normalize(self.spring_f_vec), step)
        self.vertices[0].set_coordinates(vector_add(self.vertices[0].coordinates, displacement))
        self.vertices[1].set_coordinates(vector_add(self.vertices[1].coordinates, displacement))

    # def lateral_vec(self, vector, vertex):
    #     """ returns the component of vector which is orthogonal, or lateral, to self.axis_vector """
    #     axis_vec = self.axis_vector(vertex)
    #     cos_alpha = (dot_product(vector, axis_vec) / (vector_mag(vector) * vector_mag(axis_vec)))
    #     vector_long_mag = vector_mag(vector) * cos_alpha
    #     vector_long = vector_scalar_multiply(axis_vec, vector_long_mag)
    #     vector_lat = vector_add(vector, vector_scalar_multiply(vector_long, -1))  # todo pretty up with np.array
    #     return vector_lat

    @property
    def spring_f_vec(self):
        """ returns vector sum of all spring forces """
        return vector_add(self.vertices[0].f_vector, self.vertices[1].f_vector)

    @property
    def spring_f_mag(self):
        return sum([e ** 2 for e in self.spring_f_vec]) ** 0.5


class Tendon(Member):
    """ Tension member. Exerts a pull force on its vertices equal to elongation * spring constant"""

    def __init__(self, vertices, spring_constant=1):
        Member.__init__(self, vertices)
        self.member_type = 'Tendon'
        self.nom_length = self.current_length  # set initial nom_length to current length
        self.spring_constant = spring_constant

    def set_spring_constant(self, spring_constant):
        self.spring_constant = spring_constant

    def set_nom_length(self, nom_length):
        self.nom_length = nom_length

    def set_force(self, force):
        self.nom_length = self.current_length - force / self.spring_constant

    def spring_f_vec(self, initial_vertex):
        """ returns the spring force vector that points away from initial_vertex"""
        return vector_scalar_multiply(self.axis_vector(initial_vertex=initial_vertex), self.spring_f_mag)

    def lateral_f_vec(self, vertex):
        """ return the lateral_f force that this tendon exerts on vertex. lateral_f force is the force perpendicular to
        the strut. The method is to find the projection of spring_f_vec onto strut_axis, strut_axis_f_vec, and then
        subtract strut_axis_f_vec from spring_f_vec to find lateral_f_vec"""
        strut_axis = vertex.strut.axis_vector(initial_vertex=vertex)
        strut_axis_f_vec = vector_projection(self.spring_f_vec(vertex), strut_axis)
        return np.subtract(np.array(self.spring_f_vec(vertex)), np.array(strut_axis_f_vec))

    @property
    def stretch_len(self):
        """ current_length - nom_length"""
        return self.current_length - self.nom_length

    @property
    def spring_f_mag(self):
        if self.current_length > self.nom_length:
            mag = self.spring_constant * (self.current_length - self.nom_length)
        else:
            mag = 0
        return mag

    # def unit_force_vector(self, vertex):
    #     """ return a free vector parallel to the strut that points in the direction of the force"""
    #     return normalize(self.f_vector(vertex))
    #
    def f_vector(self, vertex):
        """ return a vector parallel to the tendon that points away from vertex arg"""
        if vertex is self.vertices[1]:
            # return vector_scalar_multiply(self.vertices[0] - self.vertices[1], self.spring_f_mag)
            return vector_scalar_multiply(normalize(self.vertices[0] - self.vertices[1]), self.spring_f_mag)
        elif vertex is self.vertices[0]:
            return vector_scalar_multiply(normalize(self.vertices[1] - self.vertices[0]), self.spring_f_mag)
        else:
            raise Exception('Expected vertex to belong to member')


class Tensegrity:

    def __init__(self, coordinates, strut_vertices, tendon_vertices):
        self.vertices = [Vertex(coords) for coords in coordinates]
        self.struts = [Strut([self.vertices[vtx_indices[0]], self.vertices[vtx_indices[1]]])
                       for vtx_indices in strut_vertices]
        self.tendons = [Tendon([self.vertices[vtx_indices[0]], self.vertices[vtx_indices[1]]])
                        for vtx_indices in tendon_vertices]
        self.populate_members()

    @property
    def members(self):
        return self.struts + self.tendons

    def set_xalglib_forces(self, verbosity=0):
        """ Used at tensegrity initialization and to update after any changes to vertex coordinates"""
        xalglib_forces = self.xalglib_member_forces
        if xalglib_forces:
            if len(xalglib_forces) != len(self.members):
                raise Exception('expected xalglib_forces to be same length as self.members')
            for member, force in zip(self.members, xalglib_forces):
                member.set_xalglib_force(force)
        else:
            if verbosity > 0:
                print('****** Warning! xalglib_forces == None ******')

    def populate_members(self):
        """ Populate each vertex's list of tendons and strut (members) and store xalglib forces in each member"""
        self.set_xalglib_forces()
        for member in self.members:
            for vertex in member.vertices:
                if isinstance(member, Tendon):
                    vertex.add_tendon(member)
                elif isinstance(member, Strut):
                    vertex.set_strut(member)
                else:
                    raise Exception('expected member to be instance of Strut or Tendon')

    @property
    def vertex_array(self):
        return np.array([vertex.coordinates for vertex in self.vertices])

    @property
    def strut_array(self):
        return [[self.vertices.index(strut.vertices[0]), self.vertices.index(strut.vertices[1])]
                for strut in self.struts]

    @property
    def tendon_array(self):
        return [[self.vertices.index(tendon.vertices[0]), self.vertices.index(tendon.vertices[1])]
                for tendon in self.tendons]

    @property
    def xalglib_member_forces(self, verbosity=0):
        if verbosity > 0:
            verbose = True
        else:
            verbose = False
        tensions = o1.solve_tensegrity_tensions(self.strut_array, self.tendon_array, self.vertex_array, verbose=verbose)
        return o1.solve_tensegrity_tensions(self.strut_array, self.tendon_array, self.vertex_array, verbose=verbose)

    @property
    def xalglib_vertex_forces(self):
        """ returns the net force on each vertex based on xalglib"""
        force_vector_list = []
        i_xalglib = 0
        for vertex in self.vertices:
            force_vector = np.array([0, 0, 0])
            for member in vertex.members:
                force_vector = vector_add(force_vector, member.xalglib_force_vector(vertex))
            i_xalglib += 1
            force_vector_list.append(force_vector)
        return force_vector_list

    def set_nom_tendon_lengths(self, factor):
        """ factor is typically less than 1.0"""
        for tendon in self.tendons:
            # tendon.set_nom_length(factor)
            tendon.set_nom_length(tendon.current_length * factor)

    def equilibrium(self, err_tol=0.01):
        return False not in [strut.spring_f_mag < err_tol for strut in self.struts]

    def print_xalglib_forces(self, verbosity=0):
        forces = self.xalglib_member_forces(verbosity=verbosity)
        print('xalglib_member_forces', self.xalglib_member_forces(verbosity=verbosity))
        print('xalglib_member_forces', self.xalglib_member_forces(verbosity=2))

    def print_cylindrical(self, err_tol=0.01, vertices=True, tendons=True, struts=True, verbosity=2):
        np.set_printoptions(formatter={'float': '{: 10.3f}'.format})
        type_width = len('vertex')
        coord_width = 10  # make this equal to the width specifier in set_printoptions call above
        array_3_width = 3 * coord_width + 4  # + 4 to account for two space and two brackets
        number_width = 11
        precision = 3
        label_vertex = 'Vertex'
        label_tendon = 'Tendon'
        label_strut = 'Strut'
        label_element = 'Element'
        label_width = max([len(e) for e in [label_vertex, label_tendon, label_strut, label_element]])
        heading_coord = 'Coordinates'
        heading_v0 = f'V0 {heading_coord}'
        heading_v1 = f'V1 {heading_coord}'
        if verbosity > 0:
            print('* equilibrium is', self.equilibrium(err_tol=err_tol), 'err_tol = ', err_tol, '*')
            if vertices:
                print(f'{"Element": <{label_width}}',
                      f'{heading_coord: <{array_3_width}}',
                      f'{"Force": <{array_3_width}}')
                print(f'{" ": <{label_width}}',
                      f'{"r   ": >{coord_width}}',
                      f'{"theta": >{coord_width}}',
                      f'{"z   ": >{coord_width}}',
                      " ",
                      f'{"r   ": >{coord_width}}',
                      f'{"theta": >{coord_width}}',
                      f'{"z   ": >{coord_width}}')
                for vertex in self.vertices:
                    print(f'{label_vertex: <{label_width}}',
                          f'{str(np.array(vertex.cyl_coordinates_deg)): <{array_3_width}}',
                          f'{str(np.around(np.array(vertex.f_vector), precision)): <{array_3_width}}')
            if tendons:
                print(f'{"Element": <{label_width}}',
                      f'{heading_v0: <{array_3_width}}',
                      f'{heading_v1: <{array_3_width}}',
                      f'{"Force": <{array_3_width}}')
                print(f'{" ": <{label_width}}',
                      f'{"r   ": >{coord_width}}',
                      f'{"theta": >{coord_width}}',
                      f'{"z   ": >{coord_width}}',
                      " ",
                      f'{"r   ": >{coord_width}}',
                      f'{"theta": >{coord_width}}',
                      f'{"z   ": >{coord_width}}')
                for tendon in self.tendons:
                    vertex0 = tendon.vertices[0]
                    vertex1 = tendon.vertices[1]
                    print(f'{"Tendon": <{label_width}}',
                          f'{str(np.array(vertex0.cyl_coordinates_deg)): <{array_3_width}}',
                          f'{str(np.array(vertex1.cyl_coordinates_deg)): <{array_3_width}}',
                          f'{str(round(tendon.spring_f_mag, precision)): <{number_width}}'
                          )
            if struts:
                print(f'{"Element": <{label_width}}',
                      f'{heading_v0: <{array_3_width}}',
                      f'{heading_v1: <{array_3_width}}',
                      f'{"Twist": <{number_width}}',
                      f'{"Force Vector": <{array_3_width}}',
                      )
                print(f'{" ": <{label_width}}',
                      f'{"r   ": >{coord_width}}',
                      f'{"theta": >{coord_width}}',
                      f'{"z   ": >{coord_width}}',
                      " ",
                      f'{"r   ": >{coord_width}}',
                      f'{"theta": >{coord_width}}',
                      f'{"z   ": >{coord_width}}')
                for strut in self.struts:
                    vertex0_cyl = strut.vertices[0].cyl_coordinates_deg
                    vertex1_cyl = strut.vertices[1].cyl_coordinates_deg
                    print(f'{"Strut": <{label_width}}',
                          f'{str(np.array(vertex0_cyl)): <{array_3_width}}',
                          f'{str(np.array(vertex1_cyl)): <{array_3_width}}',
                          f'{str(round(vertex0_cyl[1] - vertex1_cyl[1], precision)): <{number_width}}',
                          f'{str(np.array(strut.spring_f_vec)): <{array_3_width}}')

    def print_spring_forces(self, err_tol, vertices=True, members=True, unit_vectors=False, verbosity=2):
        """ verbosity=0 produces no output,
            verbosity=1 produces only equilibrium status,
            verbosity 2 produces all requested information"""
        type_width = len('vertex')
        vector_width = 30
        number_width = 11
        precision = 3
        np.set_printoptions(formatter={'float': '{: 8.3f}'.format})
        label_vertex = 'Vertex'
        label_tendon = 'Tendon'
        label_strut = 'Strut'
        label_element = 'Element'
        label_width = max([len(e) for e in [label_vertex, label_tendon, label_strut, label_element]])
        heading_coord = 'Coordinates'
        heading_v0 = f'V0 {heading_coord}'
        heading_v1 = f'V1 {heading_coord}'
        heading_force = 'Force'
        if verbosity < 2:
            vertices = False
            members = False
            unit_vectors = False
        if verbosity > 0:
            if self.equilibrium(err_tol=err_tol):
                print('*** equilibrium is True ***')
            else:
                print('*** equilibrium is False ***')
            print('** spring forces **')
            if vertices:
                print(f'{"Element": <{label_width}}',
                      f'{heading_coord: <{vector_width}}',
                      f'{heading_force: <{vector_width}}')
                for vertex in self.vertices:
                    print(f'{label_vertex: <{label_width}}',
                          f'{str(np.around(np.array(vertex.coordinates), precision)): <{vector_width}}',
                          f'{str(np.around(np.array(vertex.f_vector), precision)): <{vector_width}}')
            if members:
                print(f'{"Element": <{label_width}}',
                      f'{heading_v0: <{vector_width}}',
                      f'{heading_v1: <{vector_width}}',
                      f'{"Force Mag": <{number_width}}',
                      f'{"Stretch Len": <{number_width}}'
                      )
                for tendon in self.tendons:
                    print(f'{tendon.member_type: <{label_width}}',
                          f'{str(np.around(np.array(tendon.vertices[0].coordinates), precision)): <{vector_width}}',
                          f'{str(np.around(np.array(tendon.vertices[1].coordinates), precision)): <{vector_width}}',
                          f'{str(round(tendon.spring_f_mag, precision)): <{number_width}}',
                          f'{str(round(tendon.stretch_len, precision)): <{number_width}}'
                          )
                print(f'{"Element": <{label_width}}',
                      f'{heading_v0: <{vector_width}}',
                      f'{heading_v1: <{vector_width}}',
                      f'{"Force Vec": <{number_width}}',
                      f'{"Length": <{number_width}}'
                      )
                for strut in self.struts:
                    print(f'{strut.member_type: <{label_width}}',
                          f'{str(np.around(np.array(strut.vertices[0].coordinates), precision)): <{vector_width}}',
                          f'{str(np.around(np.array(strut.vertices[1].coordinates), precision)): <{vector_width}}',
                          f'{str(np.around(np.array(strut.spring_f_vec), precision)): <{vector_width}}',
                          f'{str(round(strut.current_length, precision)): <{number_width}}'
                          )
            if unit_vectors:
                print('* Member unit vectors from simpleTools')
                for member in self.members:
                    print(f'{member.member_type: <{type_width}}',
                          f'{str(np.around(np.array(member.vertices[0].coordinates), precision)): <{vector_width}}',
                          f'{str(np.around(np.array(member.vertices[1].coordinates), precision)): <{vector_width}}',
                          # f'{str(np.around(np.array(member.unit_force_vector(member.vertices[0])), round_param)): <{vector_width}}',
                          # f'{str(np.around(np.array(member.unit_force_vector(member.vertices[0])), round_param)): <{vector_width}}'
                          )

    def solver(self, err_tol=0.01, initial_step=0.005, max_step_count=1000, verbose=False):
        """ some cases require a translation and others a rotation. All forces described here are, or result from,
        Tendon spring forces. Vertex forces are lateral_f if the angle between the vertex force and the strut
        is between pi/4 and 3pi/4. Vertex forces that lie between 0 and pi/4 or 3pi/4 and pi are longitudinal.
        If an operation results in a new vertex force that is in the opposite direction as the old vertex force
        (old dot new < 0) then an overshoot is diagnosed and the step size for that vertex is reduced"""
        step_count = 0
        # todo: add overshoot detection and step reduction so that the initial step size does not have to be smaller
        for vertex in self.vertices:
            vertex.set_step(initial_step)

        while not self.equilibrium(err_tol=err_tol) and step_count < max_step_count:
            for strut in self.struts:
                # Both vertex forces lateral_f and more parallel than orthogonal:
                # Translate strut along sum of vertex forces
                if (strut.vertices[0].f_is_lateral and strut.vertices[1].f_is_lateral and
                        abs(dot_product(normalize(strut.vertices[0].f_vector),
                                        normalize(strut.vertices[1].f_vector))) >= math.cos(math.pi / 4)):
                    strut.translate(min([vertex.step for vertex in strut.vertices]))
                    if verbose:
                        print('Both vertex forces lateral_f and more parallel than orthogonal: Translating')
                # Both vertex forces lateral_f and more orthogonal than parallel:
                # Rotate vertex with largest force about the vertex with smallest force
                elif (strut.vertices[0].f_is_lateral and strut.vertices[1].f_is_lateral and
                      abs(dot_product(normalize(strut.vertices[0].f_vector),
                                      normalize(strut.vertices[1].f_vector))) < math.cos(math.pi / 4)):
                    # vertex 0 has largest force
                    if vector_mag(strut.vertices[0].f_vector) >= vector_mag(strut.vertices[1].f_vector):
                        c_vertex = strut.vertices[1]
                        r_vertex = strut.vertices[0]
                        angle = math.asin(r_vertex.step / strut.current_length)
                        r_vertex.rotate(center=c_vertex,
                                        axis=cross_product(normalize(r_vertex.f_vector),
                                                           strut.axis_vector(r_vertex)),
                                        angle=math.asin(r_vertex.step / strut.current_length))
                    # vertex 1 has largest force
                    elif vector_mag(strut.vertices[1].f_vector) > vector_mag(strut.vertices[0].f_vector):
                        c_vertex = strut.vertices[0]
                        r_vertex = strut.vertices[1]
                        r_vertex.rotate(center=c_vertex,
                                        axis=cross_product(normalize(r_vertex.f_vector),
                                                           strut.axis_vector(r_vertex)),
                                        angle=math.asin(r_vertex.step / strut.current_length))
                    if verbose:
                        print('Both vertex forces lateral_f and more orthogonal than parallel: Rotating')
                # One vertex force is lateral_f and one is longitudinal and the lateral_f force is largest:
                # Rotate vertex with the largest force about the vertex with the smallest force
                # vertex 0 is lateral_f and largest
                elif (strut.vertices[0].f_is_lateral and strut.vertices[1].f_is_longitudinal and
                      vector_mag(strut.vertices[0].f_vector) > vector_mag(strut.vertices[1].f_vector)):
                    c_vertex = strut.vertices[1]
                    r_vertex = strut.vertices[0]
                    r_vertex.rotate(center=c_vertex,
                                    axis=cross_product(normalize(r_vertex.f_vector),
                                                       strut.axis_vector(r_vertex)),
                                    angle=math.asin(r_vertex.step / strut.current_length))
                    if verbose:
                        print('One vertex force is lateral_f and one is longitudinal and the lateral_f force',
                              'is largest: Rotating vertex 0')
                    # vertex 1 is lateral_f and largest
                elif (strut.vertices[1].f_is_lateral and strut.vertices[0].f_is_longitudinal and
                      vector_mag(strut.vertices[1].f_vector) > vector_mag(strut.vertices[0].f_vector)):
                    c_vertex = strut.vertices[0]
                    r_vertex = strut.vertices[1]
                    r_vertex.rotate(center=c_vertex,
                                    axis=cross_product(normalize(r_vertex.f_vector),
                                                       strut.axis_vector(r_vertex)),
                                    angle=math.asin(r_vertex.step / strut.current_length))
                    if verbose:
                        print('One vertex force is lateral_f and one is longitudinal and the lateral_f force ',
                              'is largest: Rotating vertex 1')
                # One vertex force is lateral_f and one is longitudinal and the longitudinal force is largest:
                # Translate along the largest force
                # vertex 0 is longitudinal and largest
                elif (strut.vertices[0].f_is_longitudinal and strut.vertices[1].f_is_lateral and
                      vector_mag(strut.vertices[0].f_vector) > vector_mag(strut.vertices[1].f_vector)):
                    strut.translate(strut.vertices[0].step)
                    if verbose:
                        print(
                            'One vertex force is lateral_f and one is longitudinal and the ',
                            'longitudinal force is largest: Translating along vertex 0 force')
                # vertex 1 is longitudinal and largest
                elif (strut.vertices[1].f_is_longitudinal and strut.vertices[0].f_is_lateral and
                      vector_mag(strut.vertices[1].f_vector) > vector_mag(strut.vertices[0].f_vector)):
                    strut.translate(strut.vertices[1].step)
                    if verbose:
                        print(
                            'One vertex force is lateral_f and one is longitudinal and the longitudinal force is largest',
                            ' Translating along vertex 0 force')
                # Both forces are longitudinal: Translate along sum of the vertex forces
                elif strut.vertices[0].f_is_longitudinal and strut.vertices[1].f_is_longitudinal:
                    strut.translate(min([vertex.step for vertex in strut.vertices]))
                    if verbose:
                        print('Both forces are longitudinal: Translating along sum of forces')
                else:
                    raise Exception('No operation found, solver cannot continue')
            step_count += 1
        if verbose:
            print('>> solver used ', step_count, ' steps')

    def plot(self, struts=True, tendons=True, lateral_f=False, vertex_f=False):
        """ plots tendons and struts using matplotlib """
        # fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.set_xlim(-5, 5)
        ax.set_ylim(-5, 5)
        ax.set_zlim(0, 10)
        if struts:
            for strut in self.struts:
                x = [vertex.coordinates[0] for vertex in strut.vertices]
                y = [vertex.coordinates[1] for vertex in strut.vertices]
                z = [vertex.coordinates[2] for vertex in strut.vertices]
                ax.plot3D(x, y, z, 'red', linewidth=3)
        if tendons:
            for tendon in self.tendons:
                x = [vertex.coordinates[0] for vertex in tendon.vertices]
                y = [vertex.coordinates[1] for vertex in tendon.vertices]
                z = [vertex.coordinates[2] for vertex in tendon.vertices]
                ax.plot3D(x, y, z, 'grey', linewidth=1)
        if lateral_f:
            # scale the vectors so the plot looks better
            scale_factor = 0.5

            for tendon in self.tendons:
                for vertex in tendon.vertices:
                    coords = [c * scale_factor for c in tendon.lateral_f_vec(vertex)]
                    ax.quiver(*(vertex.coordinates + coords), color='blue')
        if vertex_f:
            for vertex in self.vertices:
                # scale the vectors so the plot looks better
                scale_factor = 0.02
                scale_factor = 0.1
                coords = [c * scale_factor for c in vertex.f_vector]
                # ax.quiver(*(vertex.coordinates + list(vertex.f_vector)), color='green',  arrow_length_ratio=0.01)
                ax.quiver(*(vertex.coordinates + coords), color='green')
                # ax.quiver(*(vertex.coordinates + list(vertex.f_vector)), color='green')



        plt.show()
# end class Tensegrity


class TArbitrary(Tensegrity):
    """ used for testing. Allows calling function to fully describe Tensegrity """

    def __init__(self, coordinates, strut_vertices, tendon_vertices, nom_tendon_lengths):
        self.strut_vertices = strut_vertices
        self.tendon_vertices = tendon_vertices
        Tensegrity.__init__(self, coordinates, self.strut_vertices, self.tendon_vertices)
        self.set_nom_tendon_lengths(nom_tendon_lengths)


class Kite(Tensegrity):

    def __init__(self):
        coordinates = [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
        self.strut_vertices = [[0, 2], [1, 3]]
        self.tendon_vertices = [[0, 1], [1, 2], [2, 3], [3, 0]]
        Tensegrity.__init__(self, coordinates, self.strut_vertices, self.tendon_vertices)
        self.set_nom_tendon_lengths(0.5)


class PinnedKite(Tensegrity):
    """ strut are attached by a pivot at their midpoints and the angle between strut (theta) is changed to reach
    a force equilibrium """

    def __init__(self):
        coordinates = [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
        self.strut_vertices = [[0, 2], [1, 3]]
        self.tendon_vertices = [[0, 1], [1, 2], [2, 3], [3, 0]]
        Tensegrity.__init__(self, coordinates, self.strut_vertices, self.tendon_vertices)
        self.set_nom_tendon_lengths(0.9)

    def set_strut_theta(self, theta):
        """ assumes that both strut have their midpoints at [0, 0, 0] and have a length of 2"""
        self.struts[1].vertices[0].set_coordinates([math.cos(theta), math.sin(theta), 0])
        self.struts[1].vertices[1].set_coordinates([-math.cos(theta), -math.sin(theta), 0])


class Prism(Tensegrity):
    """ Intended to be used as demo with solver """

    def __init__(self, n):
        self.n = n  # strut count
        self.strut_len = 10
        self.bot_radius = 6
        self.top_radius = 6
        # the top waist polygon is twisted relative to the bottom polygon
        self.right_twist = True
        if self.right_twist:
            twist_sign = 1
        else:
            twist_sign = -1
        twist = twist_sign * (math.pi / 2 - math.pi / self.n)
        self.bot_coordinates = []
        self.top_coordinates = []
        self.s_vertices = []
        self.bot_t_vertices = []
        self.top_t_vertices = []
        self.vertical_t_vertices = []
        theta_step = 2 * math.pi / self.n
        theta = 0
        z_bot = 0
        z_top = self.strut_len
        for i in range(self.n):
            self.bot_coordinates.append(cartesian_coordinates([self.bot_radius, theta, z_bot]))
            self.top_coordinates.append(cartesian_coordinates([self.top_radius, theta + twist, z_top]))
            self.bot_t_vertices.append([i, (i + 1) % self.n])
            self.top_t_vertices.append([i + n, (i + 1) % self.n + self.n])
            self.vertical_t_vertices.append([i, i + self.n])
            self.s_vertices.append([i, (i + twist_sign) % self.n + self.n])
            theta += theta_step
        Tensegrity.__init__(self, self.bot_coordinates + self.top_coordinates, self.s_vertices,
                            self.bot_t_vertices + self.top_t_vertices + self.vertical_t_vertices)
        self.bot_vertices = self.vertices[0:len(self.bot_coordinates)]
        self.top_vertices = self.vertices[len(self.bot_coordinates):]
        self.bot_tendons = self.tendons[0:len(self.bot_t_vertices)]
        self.top_tendons = self.tendons[len(self.bot_t_vertices):len(self.top_t_vertices) + len(self.bot_t_vertices)]
        self.vertical_tendons = self.tendons[len(self.bot_t_vertices) + len(self.top_t_vertices):]

    def balance_forces(self, verbose=0, debug=True):
        """ Assumes a known good symmetrical prism.
        Methodology:
        1. Assign an arbitrary force equally to bottom waist tendons
        2. Find the vertical tendon force that will balance the waist tendons in the plane orthogonal to the strut and
        containing the vertex
        3. Find the top waist forces that will balance the vertical tendon forces in the top orthogonal planes
        4. Check for equilibrium
        """
        # 1. Assign an arbitrary force equally to bottom waist tendons
        bot_t_force = 10
        for tendon in self.bot_tendons:
            tendon.set_force(bot_t_force)
        # 2. Find the vertical tendon force that will balance the waist tendons in the plane orthogonal to the strut and
        # containing the vertex
        for vertex in self.bot_vertices:
            waist_lateral_f_vecs = []
            for tendon in vertex.tendons:
                if tendon in self.bot_tendons:
                    waist_lateral_f_vecs.append(tendon.lateral_f_vec(vertex))
            if len(waist_lateral_f_vecs) != 2:
                raise Exception('Expected to find 2 bot_tendons for bot_vertex')
            waist_lateral_f_vec = vector_add(waist_lateral_f_vecs[0], waist_lateral_f_vecs[1])
            vertical_lateral_f_vec = vector_scalar_multiply(waist_lateral_f_vec, -1)
            vertical_tendon = [tendon for tendon in vertex.tendons if tendon in self.vertical_tendons][0]
            if vector_mag(cross_product(waist_lateral_f_vec, vertical_tendon.lateral_f_vec(vertex))) > stable_tol:
                raise Exception('This tensegrity is not stable')
            alpha = vec_angle(vertical_tendon.axis_vector(vertex), vertex.strut.axis_vector(vertex))
            # phi = vec_angle(vertical_tendon.axis_vector(vertex), waist_lateral_f_vec)
            phi = vec_angle(vertical_tendon.axis_vector(vertex), vertical_lateral_f_vec)
            # cos_alpha = (dot_product(vertical_tendon.axis_vector(vertex), vertex.strut.axis_vector(vertex)) /
            #              (vector_mag(vertical_tendon.axis_vector(vertex)) *
            #               vector_mag(vertex.strut.axis_vector(vertex))))
            # sin_alpha = (1 - cos_alpha ** 2) ** 0.5
            vertical_tendon.set_force(vector_mag(vertical_lateral_f_vec) / math.cos(phi))
            # vertical_tendon.set_force(vector_mag(waist_lateral_f_vec) / abs(vec_cos(vertical_tendon.axis_vector(vertex),
            #                                                                     waist_lateral_f_vec)))

            # vertical_tendon.set_force(vector_mag(waist_lateral_f_vec) * math.sin(alpha))

            if debug:
                pass
                # print('dot_product(waist_lateral_f_vecs[0])',
                #       dot_product(waist_lateral_f_vecs[0], vertex.strut.axis_vector(vertex)))
                # print('dot_product(waist_lateral_f_vecs[1])',
                #       dot_product(waist_lateral_f_vecs[1], vertex.strut.axis_vector(vertex)))
                # print('dot_product(waist_lateral_f_vec)',
                #       dot_product(waist_lateral_f_vec, vertex.strut.axis_vector(vertex)))
                # print('** cross waist_lateral_f_vec x vertical_lateral_f_vec',
                #       cross_product(waist_lateral_f_vec, vertical_tendon.lateral_f_vec(vertex)))
                # print('waist_lateral_f_vecs[0] mag', vector_mag(waist_lateral_f_vecs[0]))
                # print('waist_lateral_f_vecs[1] mag', vector_mag(waist_lateral_f_vecs[1]))
                # waist_sum = vector_add(*[tendon.lateral_f_vec(vertex) for
                #                          tendon in vertex.tendons if tendon in self.bot_tendons])
                # print('lateral vec sum', vector_add(vertical_tendon.lateral_f_vec(vertex), waist_sum))
            if verbose > 1:
                print('bot vertices f_vector', vertex.f_vector, 'cross with strut axis',
                      cross_product(vertex.f_vector, vertex.strut.axis_vector(vertex)))
        # 3. Find the top waist forces that will balance the vertical tendon forces in the top orthogonal planes
        for vertex in self.top_vertices:
            vertical_tendon = [tendon for tendon in vertex.tendons if tendon in self.vertical_tendons][0]
            # debug
            lat_f_vec = vertical_tendon.lateral_f_vec(vertex)
            waist_lat_sum_vec = vector_scalar_multiply(vertical_tendon.lateral_f_vec(vertex), -1)
            alpha = vec_angle(*[tendon.axis_vector(vertex) for tendon in vertex.tendons if tendon in self.top_tendons])
            waist_lat_f_mag = (vector_mag(waist_lat_sum_vec) / 2) / math.cos(alpha / 2)
            for tendon in vertex.tendons:
                if tendon in self.top_tendons:
                    theta = vec_angle(tendon.axis_vector(vertex), vertex.strut.axis_vector(vertex))
                    tendon.set_force(waist_lat_f_mag / math.sin(theta))
            if verbose > 1:
                print('top vertices f_vector', vertex.f_vector, 'cross with strut axis',
                      cross_product(vertex.f_vector, vertex.strut.axis_vector(vertex)))
        # 4. Check for equilibrium
        if verbose > 1:
            print('equilibrium is', self.equilibrium(0.5))

    # def solver_strut_twister(self):
    #     """ attempts to reach equilibrium by changing position of the top of the struts. Keeps all radial symmetries"""
    #     # todo: take into account the fact that err_tol is Newtons and step size is in radians
    #     # perhaps use feedback loop to modify step size in order to get the desired Newton step size?
    #     err_tol = 0.1
    #     max_steps = 1000
    #     # step_sizes = [5 * err_tol, 2 * err_tol, 0.9 * err_tol]
    #     phi_step_size = 1.5 * err_tol  # wild guess
    #     theta_step_size = phi_step_size
    #     phi_direction = -1
    #     theta_direction = 1
    #     # i_step_size = 0
    #     step_count = 0
    #     while not self.equilibrium(err_tol=err_tol):
    #         # modify phi (spherical polar angle) of top vertex in all struts equally
    #         for strut in self.struts:
    #             # strut.modify_phi(phi_direction * step_sizes[i_step_size])
    #             initial_strut_force_mag = self.struts[0].spring_force_magnitude
    #             strut.modify_phi(phi_direction * step_size)
    #             if self.struts[0].spring_force_magnitude > initial_strut_force_mag:
    #                 # if we are going the wrong way then reverse direction
    #                 phi_direction = -phi_direction
    #                 # phi_better = False
    #             elif self.struts[0].spring_force_magnitude == initial_strut_force_mag:
    #                 raise Exception('no change in spring force found after change in phi')
    #             initial_strut_force_mag = self.struts[0].spring_force_magnitude
    #             strut.modify_theta(theta_direction * step_size)
    #
    #             # else:
    #             # phi_better = True
    #         step_count += 1
    #         if step_count > max_steps:
    #             raise Exception('max step count exceeded')
    #         # modify theta (spherical azimuthal angle) of top vertex in all struts equally


class PrismNeq3(Tensegrity):
    """ N = 3 prism from olof_1.py"""

    def __init__(self):
        alpha = 5 * math.pi / 6
        self.coordinates = [[math.cos(0), math.sin(0), 0],
                            [math.cos(2 * math.pi / 3), math.sin(2 * math.pi / 3), 0],
                            [math.cos(4 * math.pi / 3), math.sin(4 * math.pi / 3), 0],
                            [math.cos(alpha), math.sin(alpha), 1],
                            [math.cos(alpha + 2 * math.pi / 3), math.sin(alpha + 2 * math.pi / 3), 1],
                            [math.cos(alpha + 4 * math.pi / 3), math.sin(alpha + 4 * math.pi / 3), 1]]
        # these numbers are indices in the vertices array for the vertices on either end of the member
        self.strut_vertices = [[0, 3],
                               [1, 4],
                               [2, 5]]
        self.tendon_vertices = [[0, 1],
                                [1, 2],
                                [2, 0],
                                [3, 4],
                                [4, 5],
                                [5, 3],
                                [1, 3],
                                [2, 4],
                                [0, 5]]
        Tensegrity.__init__(self, self.coordinates, self.strut_vertices, self.tendon_vertices)

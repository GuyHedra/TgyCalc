import math
import olof_1 as o1
from numpy import array, around
import numbers
from scipy.spatial.transform import Rotation as R

""" simpleTools.py is a simplified version of polyTools intended for use within TgyCalc"""


def vector_add(v0, v1):
    return [c0 + c1 for c0, c1 in zip(v0, v1)]


def vector_mag(vector):
    """ return magnitude of a vector"""
    return sum([c ** 2 for c in vector]) ** 0.5


def vector_scalar_multiply(vector, scalar):
    if isinstance(vector, list) and isinstance(scalar, numbers.Number):
        return [element * scalar for element in vector]
    else:
        raise Exception('vector must be a list and scalar must be a number')


def dot_product(v0, v1):
    """returns the dot product of two vectors"""
    return sum([v0_coord * v1_coord for v0_coord, v1_coord in zip(v0, v1)])


def cross_product(vector0, vector1):
    """Returns the cross product between two free vectors"""
    return [vector0[1] * vector1[2] - vector0[2] * vector1[1],
            vector0[2] * vector1[0] - vector0[0] * vector1[2],
            vector0[0] * vector1[1] - vector0[1] * vector1[0]]


def normalize_vector(vector):
    magnitude = sum([element ** 2 for element in vector]) ** 0.5
    if magnitude != 0:
        return [element / magnitude for element in vector]
    else:
        return vector


def cartesian_coordinates(cyl_coords):
    """ cylindrical are [rho(r), phi(azimuth), z]"""
    rho, phi, z = cyl_coords
    return [rho * math.cos(phi), rho * math.sin(phi), z]


# def axis_angle_to_quaternion(axis, theta):
#     """ Convert from axis-angle rotations to quaternions."""
#     x, y, z = normalize_vector(axis)
#     theta /= 2
#     w = math.cos(theta)
#     x = x * math.sin(theta)
#     y = y * math.sin(theta)
#     z = z * math.sin(theta)
#     return Quaternion((w, x, y, z))


# class Quaternion:
#     """ Quaternion is used to store and operate on quaternions """
#     def __init__(self, coordinates):
#         self.coordinates = coordinates
#         # tuple.__init__(self)
#         # if (isinstance(quaternion, tuple) and len(quaternion) == 4
#         #         and False not in [isinstance(coordinate, numbers.Number) for coordinate in quaternion]):
#         #     self.coords = quaternion
#         # else:
#         #     raise TypeError('quaternion argument must be a 4 element tuple of Number')
#
#     @property
#     def angle(self):
#         return self.coordinates[0]
#
#     @property
#     def vector(self):
#         return Vector(self.coordinates[1:])
#
#     @property
#     def conjugate(self):
#         """ Return the conjugate of a quaternion"""
#         # w, x, y, z = self
#         return Quaternion([self.angle, -self.coordinates[0], -self.coordinates[1], -self.coordinates[2]])
#
#     def __eq__(self, other):
#         if not isinstance(other, tuple) and len(other) == 4:
#             raise TypeError('bad operand for Quaternion.__eq__', other.__class__)
#         if False not in [abs(s_coord - o_coord) < tolerance
#                          for s_coord, o_coord in zip(self.coordinates, other)]:
#             return True
#         else:
#             return False
#
#     def __mul__(self, other):
#         if isinstance(other, Quaternion):
#             # Return the product of two quaternions
#             w1, x1, y1, z1 = self.coordinates
#             w2, x2, y2, z2 = other.coordinates
#             w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
#             x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
#             y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
#             z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
#             return Quaternion((w, x, y, z))
#         # elif isinstance(other, Point):
#         #     # return the product of this quaternion and a free vector (rotate the free vector by the quaternion)
#         #     q = Quaternion((0.0,) + other.coords)
#         #     return Point((self * q * self.conjugate)[1:])
#         elif isinstance(other, FreeVector):
#             # return the product of this quaternion and a free vector
#             q = Quaternion([0.0] + other)
#             interim_value = self * q
#             return [interim_value * self.conjugate][1:]
#         else:
#             raise TypeError('bad argument, expected Number, Point or Quaternion')
#

# class FreeVector(list):
#     """ Vector is used to store and operate on vectors in 3-space """
#     def __init__(self, coordinates):
#         list.__init__(self)
#         """ initial_point and terminal_point must be Point objects"""
#         self.coordinates = coordinates
#
#     @property
#     def magnitude(self):
#         return sum([c ** 2 for c in self.coordinates]) ** 0.5
#
#     # @property
#     # def midpoint(self):
#     #     """ Returns the midpoint of Vector self"""
#     #     return Point(((self.points[0].coords[0] + self.points[1].coords[0]) / 2,
#     #                   (self.points[0].coords[1] + self.points[1].coords[1]) / 2,
#     #                   (self.points[0].coords[2] + self.points[1].coords[2]) / 2
#     #                   ))
#
#     # @property
#     # def free(self):
#     #     """ Return the free vector, which is this vector translated such that the initial point lies on the origin"""
#     #     return Point(self.points[1] - self.points[0])
#
#     # @property
#     # def reverse(self):
#     #     return Vector(self.points[1], self.points[0])
#
#     @property
#     def normalize(self, error_tol=1e-12):
#         """ Return a unit vector parallel to self"""
#         deltas = [t_coordinate - i_coordinate
#                   for t_coordinate, i_coordinate in zip(self.points[1], self.points[0])]
#         mag2 = sum([delta ** 2 for delta in deltas])
#         mag = mag2 ** 0.5
#         if mag < error_tol:
#             raise ValueError('magnitude of vector is less than error_tol')
#         else:
#             return Vector(self.points[0], Point(tuple([delta / mag for delta in deltas])) + self.points[0])
#
#     def __eq__(self, other):
#         """ Returns true if both initial points and terminal points are equal. Note that this definition of equal
#         is different than the mathematical definition of vector equality in which only magnitude and direction need
#         agree """
#         if not isinstance(other, Vector):
#             raise TypeError('Argument other must be of type Vector')
#         if self.points[0] == other.points[0] and self.points[1] == other.points[1]:
#             return True
#         else:
#             return False
#
#     def __add__(self, other):
#         """ adds Point other to self by translating other such that others initial point is at self's terminal point,
#         and returning a vector with an initial point equal to self's initial point and a terminal point equal to the
#         post-translation terminal point of other"""
#         if not isinstance(other, Vector):
#             raise TypeError('Argument other must be of type Vector')
#         offset = self.points[1] - other.points[0]
#         return Vector(self.points[0], other.points[1] + offset)
#
#     def __mul__(self, other):
#         """ multiply self by the scalar other"""
#         if not isinstance(other, numbers.Number):
#             raise TypeError('Argument other must be of type Number')
#         deltas = Point(tuple([(t_coordinate - i_coordinate) * other
#                               for t_coordinate, i_coordinate in zip(self.points[1], self.points[0])]))
#         return Vector(self.points[0], self.points[0] + deltas)
#
#     def __repr__(self):
#         """ returns a string representation of self"""
#         return '[' + str(self.points[0]) + ', ' + str(self.points[1]) + ']'
#
#     def angle(self, vector):
#         if not isinstance(vector, Vector):
#             raise TypeError('Argument normal_vector must be of type Point')
#         value = math.acos(dotproduct(self.free, vector.free)/(self.magnitude * vector.magnitude))
#         # if the angle is greater than pi / 2 then it is the obtuse angle, and we want the acute angle
#         if value > math.pi / 2:
#             value = math.pi - value
#         return value
#
#     def translate(self, point):
#         """ translate the Vector such that the initial point is located at the point argument"""
#         if not isinstance(point, Point):
#             raise TypeError('point argument must be an object of type Point')
#         return Vector(point, point + self.points[1] - self.points[0])
#
#     def rotate(self, axis, angle):
#         q = axis_angle_to_quaternion(axis, angle)
#         return q * self

    # def quaternion(self, initial_vector=Point((0, 0, 1)), initial_theta_vector=Point((0, 1, 0)),
    #                theta_point=Point((0, 0, 0)), debug=False, verbose=False):
    #     """
    #     Designed for use by Tendon and Strut (vector) objects, which inherit this class
    #     The initial orientation of the vector object is parallel to the z axis
    #     The final orientation of the vector object is from self.points[0] to self.points[1]
    #     The theta_point specifies the rotation about the vector self such that the initial_theta_vector
    #     is lying in the plane formed by the vector self and theta_point.
    #     Put another way, keep the side of the vector object that should be facing outward facing outward.
    #     """
    #     quat0 = rotation_quaternion(initial_vector, self, verbose=False)
    #     intermediate_theta_vector = quat0 * initial_theta_vector
    #     # intermediate_theta_vector needs to be rotated about the axis of self such that the rotated
    #     # intermediate_theta_vector lies in the plane formed by the vector self and theta_point
    #     # We define a point P lying on the line defined by the points of self. P is the closest point on the line
    #     # self to theta_point. We define a vector final_theta_vector from P to theta_point.
    #     # final_theta_vector is the vector to which we wish to rotate intermediate_theta_vector
    #     p = closest_point(self, theta_point)  # should equal midpoint for polyhedra
    #     final_theta_vector = Vector(theta_point, p)
    #     rotangle1 = rotangle(intermediate_theta_vector, final_theta_vector.free, self.free.normalize)
    #     quat1 = axisangle_to_q(self.free, rotangle1)
    #     if verbose:
    #         print('intermediate_theta_vector', intermediate_theta_vector, '\nP', p, 'midpoint', self.midpoint,
    #               '\nfinal_theta_vector', final_theta_vector, 'quat1', quat1)
    #     quat_final = quat1 * quat0
    #     if debug:
    #         print('quat_final', quat_final)
    #         self.intermediate_theta_vector = intermediate_theta_vector
    #         self.final_theta_vector = final_theta_vector
    #     return quat_final


class Vertex:
    def __init__(self, cartesian=None, cylindrical=None):
        if cylindrical:
            self.coordinates = cartesian_coordinates(cylindrical)
        elif cartesian:
            self.coordinates = cartesian
        else:
            raise Exception('Either cartesian or cylindrical must be passed')
        self.strut = None
        self.tendons = []
        self.anchor = False

    def set_coordinates(self, coordinates, cyl_coords=None):
        if cyl_coords:
            self.coordinates = cartesian_coordinates(cyl_coords)
        else:
            self.coordinates = coordinates

    def set_strut(self, strut):
        self.strut = strut

    def set_anchor(self, anchor):
        if isinstance(anchor, bool):
            self.anchor = anchor
        else:
            raise Exception('anchor must be True or False')

    def add_tendon(self, tendon):
        self.tendons.append(tendon)

    def rotate(self, center, axis, angle):
        """ Rotates self about center and axis by angle radians.
        Translates self such that center is at the origin, rotates self by angle radians about axis and then
        translates self such that center is back at its original location"""
        point = self - center
        r = R.from_rotvec(angle * array(axis))
        interim_coordinates = r.apply(point)
        self.coordinates = interim_coordinates + center
        # quat = axis_angle_to_quaternion(axis, angle)
        # self.coordinates = quat * point + center

    @property
    def members(self):
        return [self.strut] + self.tendons

    @property
    def spring_force_vector(self):
        """ returns the vector sum af all tendon spring forces acting on this vertex """
        force_vector = [0, 0, 0]
        for tendon in self.tendons:
            force_vector = vector_add(tendon.spring_force_vector(self), force_vector)
        return force_vector

    # @property
    # def total_spring_force_vector(self):
    #     #debug
    #     sfv = self.strut.spring_force_vector(self)
    #     return vector_add(self.tendon_spring_force_vector, self.strut.spring_force_vector(self))

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
        self.force = None
        self.xalglib_force = None
        self.nom_length = self.current_length  # set initial nom_length to current length
        self.spring_constant = 1

    def xalglib_force_vector(self, vertex):
        force_vector = self.unit_force_vector(vertex)  # debug
        return [element * self.xalglib_force for element in self.unit_force_vector(vertex)]

    def set_xalglib_force(self, force):
        self.xalglib_force = force

    def set_nom_length(self, nom_length):
        self.nom_length = nom_length

    @property
    def current_length(self):
        return self.vertices[0].distance(self.vertices[1])

    @property
    def anchor(self):
        return True in [vertex.anchor for vertex in self.vertices]

    @property
    def midpoint(self):
        """ Returns the midpoint of self"""
        return [(c0 + c1) / 2 for c0, c1 in zip(self.vertices[0].coordinates, self.vertices[1].coordinates)]

    def rotate(self, center_vertex, axis, angle):
        for vertex in self.vertices:
            vertex.rotate(center=center_vertex.coordinates, axis=axis, angle=angle)

    def unit_force_vector(self, vertex):
        if ((vertex is self.vertices[0] and isinstance(self, Strut)) or
                (vertex is self.vertices[1] and isinstance(self, Tendon))):
            unit_vector = normalize_vector(self.vertices[0] - self.vertices[1])
        elif ((vertex is self.vertices[1] and isinstance(self, Strut)) or
                (vertex is self.vertices[0] and isinstance(self, Tendon))):
            unit_vector = normalize_vector(self.vertices[1] - self.vertices[0])
        else:
            raise Exception('Expected vertex to belong to member')
        return unit_vector

    def spring_force_vector(self, vertex):
        # debug
        # ufv = self.unit_force_vector(vertex)
        # sfm = self.spring_force_magnitude
        # print('ufv, sfm', ufv, sfm)
        return vector_scalar_multiply(self.unit_force_vector(vertex), self.spring_force_magnitude)


class Strut(Member):

    def __init__(self, vertices):
        Member.__init__(self, vertices)
        self.member_type = 'Strut '

    def translate(self, step):
        """ sum vertex sprint vectors and translate strut by step in said vector direction"""
        displacement = vector_scalar_multiply(normalize_vector(self.spring_force_vector_sum), step)
        self.vertices[0].set_coordinates(vector_add(self.vertices[0].coordinates, displacement))
        self.vertices[1].set_coordinates(vector_add(self.vertices[1].coordinates, displacement))

    @property
    def spring_force_vector_sum(self):
        """ returns spring force if not an anchor, else zero vector"""
        return vector_add(self.vertices[0].spring_force_vector, self.vertices[1].spring_force_vector)

    @property
    def spring_force_magnitude(self):
        return sum([e ** 2 for e in self.spring_force_vector_sum]) ** 0.5

    def vector(self, vertex):
        """ return a free vector parallel to the strut that terminates at the vertex arg"""
        # todo: move this func from strut to member
        if self.vertices[0] is vertex:
            return self.vertices[0] - self.vertices[1]
        elif self.vertices[1] is vertex:
            return self.vertices[1] - self.vertices[0]
        else:
            raise Exception('strut does not contain vertex')


class Tendon(Member):
    def __init__(self, vertices):
        Member.__init__(self, vertices)
        self.member_type = 'Tendon'

    @property
    def spring_force_magnitude(self):
        if self.current_length > self.nom_length:
            mag = self.spring_constant * (self.current_length - self.nom_length)
        else:
            mag = 0
        return mag


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

    def set_xalglib_forces(self):
        """ used at tensegrity initialization and to update after any changes to vertex coordinates"""
        xalglib_forces = self.xalglib_member_forces
        if xalglib_forces:
            if len(xalglib_forces) != len(self.members):
                raise Exception('expected xalglib_forces to be same length as self.members')
            for member, force in zip(self.members, xalglib_forces):
                member.set_xalglib_force(force)
        else:
            print('****** Warning! xalglib_forces == None ******')

    def populate_members(self):
        """ populate each vertex's list of tendons and strut (members) and store xalglib forces in each member"""
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
        return array([vertex.coordinates for vertex in self.vertices])

    @property
    def strut_array(self):
        return [[self.vertices.index(strut.vertices[0]), self.vertices.index(strut.vertices[1])]
                for strut in self.struts]

    @property
    def tendon_array(self):
        return [[self.vertices.index(tendon.vertices[0]), self.vertices.index(tendon.vertices[1])]
                for tendon in self.tendons]

    @property
    def xalglib_member_forces(self):
        return o1.solve_tensegrity_tensions(self.strut_array, self.tendon_array, self.vertex_array, verbose=False)

    @property
    def xalglib_vertex_forces(self):
        """ returns the net force on each vertex based on xalglib"""
        force_vector_list = []
        # xalglib_forces = self.xalglib_member_forces
        i_xalglib = 0
        for vertex in self.vertices:
            force_vector = array([0, 0, 0])
            for member in vertex.members:
                force_vector = vector_add(force_vector, member.xalglib_force_vector(vertex))
            i_xalglib += 1
            force_vector_list.append(force_vector)
        return force_vector_list

    def set_nom_tendon_lengths(self, factor):
        """ factor is typically less than 1.0"""
        for tendon in self.tendons:
            tendon.set_nom_length(factor)

    def equilibrium(self, err_tol=0.01):
        strut_forces = [strut.spring_force_magnitude for strut in self.struts]
        eq_matrix = [strut.spring_force_magnitude < err_tol for strut in self.struts]
        ret_val = False not in [strut.spring_force_magnitude < err_tol for strut in self.struts]
        # print('>>>>>> equilibrium debug', strut_forces, eq_matrix, ret_val)
        # return True not in [strut.spring_force_magnitude > err_tol for strut in self.struts]
        return False not in [strut.spring_force_magnitude < err_tol for strut in self.struts]

    def print_spring_forces(self, err_tol, vertices=True, members=True, unit_vectors=False):
        type_width = len('vertex')
        coord_width = 20
        round_param = 3
        if self.equilibrium(err_tol=err_tol):
            print('*** This tensegrity is in equilibrium ***')
        else:
            print('*** This tensegrity is NOT in equilibrium ***')
        print('** spring forces **')
        if vertices:
            print('*      Coordinates          Tendon Force Vector')
            for vertex in self.vertices:
                print('Vertex',
                      f'{str(around(array(vertex.coordinates), round_param)): <{coord_width}}',
                      f'{str(around(array(vertex.spring_force_vector), round_param)): <{coord_width}}')
        if members:
            print('*      V0 Coordinates       V1 Coordinates       Force      Unit Vector 0        Unit Vector 1'
                  '   current length')
            for tendon in self.tendons:
                # print(f'{0: <{type_width}}'.format(tendon.member_type),
                print(f'{tendon.member_type: <{type_width}}',
                      f'{str(around(array(tendon.vertices[0].coordinates), round_param)): <{coord_width}}',
                      f'{str(around(array(tendon.vertices[1].coordinates), round_param)): <{coord_width}}',
                      round(tendon.spring_force_magnitude, round_param), '     ',
                      f'{str(around(array(tendon.unit_force_vector(tendon.vertices[0])), round_param)): <{coord_width}}',
                      f'{str(around(array(tendon.unit_force_vector(tendon.vertices[1])), round_param)): <{coord_width}}',
                      f'{str(round(tendon.current_length, round_param)): <{coord_width}}'
                      )
            for strut in self.struts:
                print(f'{strut.member_type: <{type_width}}',
                      f'{str(around(array(strut.vertices[0].coordinates), round_param)): <{coord_width}}',
                      f'{str(around(array(strut.vertices[1].coordinates), round_param)): <{coord_width}}',
                      f'{str(around(array(strut.spring_force_vector_sum), round_param)): <{coord_width}}',
                      f'{str(round(strut.current_length, round_param)): <{coord_width}}'
                      )
        if unit_vectors:
            print('* Member unit vectors from simpleTools')
            for member in self.members:
                print(f'{member.member_type: <{type_width}}',
                      f'{str(around(array(member.vertices[0].coordinates), round_param)): <{coord_width}}',
                      f'{str(around(array(member.vertices[1].coordinates), round_param)): <{coord_width}}',
                      f'{str(around(array(member.unit_force_vector(member.vertices[0])), round_param)): <{coord_width}}',
                      f'{str(around(array(member.unit_force_vector(member.vertices[0])), round_param)): <{coord_width}}'
                      )

    def solver_spring_forces(self, err_tol=0.01, initial_step=0.005, max_step_count=1000, verbose=False):
        """ If the dot product of the vertex forces is above a threshold (e.g. 0) then the forces are mostly
        translational and we will translate the strut.
        If the dot product of the vertex forces is below a threshold (e.g. 0) then the forces are mostly
        rotational and we will rotate the strut"""
        step_count = 0
        translate_threshold = 0
        # step = initial_step
        for strut in self.struts:
            strut.step = initial_step
        while not self.equilibrium(err_tol=err_tol) and step_count < max_step_count:
            for strut in self.struts:
                strut.previous_force_vector_sum = strut.spring_force_vector_sum
                # dot product less than threshold, or strut anchored, so rotate
                if (dot_product(strut.vertices[0].spring_force_vector, strut.vertices[1].spring_force_vector) <
                        translate_threshold or strut.anchor):
                    angle = math.asin(strut.step / strut.current_length)
                    if (vector_mag(strut.vertices[0].spring_force_vector) >
                            vector_mag(strut.vertices[1].spring_force_vector)) or strut.vertices[1].anchor:
                        rotation_vertex = strut.vertices[0]
                        center_vertex = strut.vertices[1]
                    else:
                        rotation_vertex = strut.vertices[1]
                        center_vertex = strut.vertices[0]
                    strut_vector = strut.vector(vertex=rotation_vertex)
                    axis = normalize_vector(cross_product(strut_vector, rotation_vertex.spring_force_vector))
                    if verbose:
                        print('>> solver rotation: center vtx:', center_vertex.coordinates, 'axis:', axis,
                              'angle', angle)
                    strut.rotate(center_vertex, axis, angle)
                else:
                    # dot product greater than threshold, so translate
                    if verbose:
                        print('>> solver translation', strut.step)
                    strut.translate(strut.step)
                if dot_product(strut.previous_force_vector_sum, strut.spring_force_vector_sum) < 0:
                    strut.step = strut.step/2
                    print('*** detected an overshoot ***')
            if verbose:
                print('step count:', step_count)
            step_count += 1
            #debug
            if step_count == 14:
                print('*****step count 14')
            if verbose:
                print('strut[0] force mag', self.struts[0].spring_force_magnitude,
                      'strut[1] force mag', self.struts[1].spring_force_magnitude)
                print('>> solver strut[0]',
                      self.struts[0].vertices[0].coordinates,
                      self.struts[0].vertices[1].coordinates,
                      'spring force sum', self.struts[0].spring_force_vector_sum
                      'strut[1]',
                      self.struts[1].vertices[0].coordinates,
                      self.struts[1].vertices[1].coordinates)
        if verbose:
            print('>> solver used ', step_count, ' steps')
# end class Tensegrity


class Kite(Tensegrity):

    def __init__(self):
        coordinates = [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
        self.strut_vertices = [[0, 2], [1, 3]]
        self.tendon_vertices = [[0, 1], [1, 2], [2, 3], [3, 0]]
        Tensegrity.__init__(self, coordinates, self.strut_vertices, self.tendon_vertices)
        self.set_nom_tendon_lengths(0.5)

    def solver_strut_0x(self):
        """ Moves struts[0] along x axis until kite is in equilibrium. Assumes the required symmetries exist """
        initial_step = 0.01
        max_step_count = 100
        step_count = 0
        while not self.equilibrium(err_tol=0.1):
            if step_count > max_step_count:
                print('**** max_step_count reached ****')
                break
            error_signal = dot_product([1, 0, 0], self.struts[0].spring_force_vector_sum)
            if error_signal < 0:
                step = -initial_step
            else:
                step = initial_step
            for vertex in self.struts[0].vertices:
                vertex.coordinates = vector_add([step, 0, 0], vertex.coordinates)
            # print('solver: spring_force_vector_sum',  self.struts[0].spring_force_vector_sum)
            step_count += 1


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
        self.bot_radius = 5
        self.top_radius = 5
        self.strut_rho = math.pi / 2
        bot_coordinates = []
        top_coordinates = []
        s_vertices = []
        bot_t_vertices = []
        top_t_vertices = []
        vertical_t_vertices = []
        # we will initialize the prism with vertical struts and let the solver find the correct strut_rho, where
        # strut_rho is the angle between the xy plane and each strut
        # the solver will use vary strut_rho until a stable solution is found
        theta_step = 2 * math.pi / self.n
        phi = 0
        z_bot = 0
        z_top = self.strut_len
        for i in range(self.n):
            bot_coordinates.append(cartesian_coordinates([self.bot_radius, phi, z_bot]))
            top_coordinates.append(cartesian_coordinates([self.bot_radius, phi, z_top]))
            bot_t_vertices.append([i, (i + 1) % self.n])
            top_t_vertices.append([i + n, (i + 1) % self.n + self.n])
            vertical_t_vertices.append([i, (i + 1) % self.n + self.n])
            s_vertices.append([i, i + self.n])
            phi += theta_step
        Tensegrity.__init__(self, bot_coordinates + top_coordinates, s_vertices,
                            bot_t_vertices + top_t_vertices + vertical_t_vertices)

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
                            [math.cos(2*math.pi/3), math.sin(2*math.pi/3), 0],
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
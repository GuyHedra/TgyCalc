import math
import olof_1 as o1
from numpy import array, around
import numbers

""" simpleTools.py is a simplified version of polyTools intended for use within TgyCalc"""


def vector_add(v0, v1):
    return [c0 + c1 for c0, c1 in zip(v0, v1)]


def vector_scalar_multiply(vector, scalar):
    if isinstance(vector, list) and isinstance(scalar, numbers.Number):
        return [element * scalar for element in vector]
    else:
        raise Exception('vector must be a list and scalar must be a number')


def dot_product(v0, v1):
    """returns the dot product of two vectors"""
    return sum([v0_coord * v1_coord for v0_coord, v1_coord in zip(v0, v1)])


def normalize_vector(vector):
    magnitude = sum([element ** 2 for element in vector]) ** 0.5
    return [element / magnitude for element in vector]

def cartesian_coordinates(cyl_coords):
    """ cylindrical are [rho(r), phi(azimuth), z]"""
    rho, phi, z = cyl_coords
    return [rho * math.cos(phi), rho * math.sin(phi), z]


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

    def set_coordinates(self, coordinates, cyl_coords=None):
        if cyl_coords:
            self.coordinates = cartesian_coordinates(cyl_coords)
        else:
            self.coordinates = coordinates

    def set_strut(self, strut):
        self.strut = strut

    def add_tendon(self, tendon):
        self.tendons.append(tendon)

    @property
    def members(self):
        return [self.strut] + self.tendons

    @property
    def tendon_spring_force_vector(self):
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
        return [self_coord - other_coord for self_coord, other_coord in zip(self.coordinates, other.coordinates)]


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

    @property
    def spring_force_vector_sum(self):
        return vector_add(self.vertices[0].tendon_spring_force_vector, self.vertices[1].tendon_spring_force_vector)

    @property
    def spring_force_magnitude(self):
        return sum([e ** 2 for e in self.spring_force_vector_sum]) ** 0.5


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

    def populate_members(self):
        """ populate each vertex's list of tendons and strut (members) and store xalglib forces in each member"""
        xalglib_forces = self.xalglib_member_forces
        if xalglib_forces:
            if len(xalglib_forces) != len(self.members):
                raise Exception('expected xalglib_forces to be same length as self.members')
            for member, force in zip(self.members, xalglib_forces):
                member.set_xalglib_force(force)
        else:
            for member in self.members:
                for vertex in member.vertices:
                    if isinstance(member, Tendon):
                        vertex.add_tendon(member)
                    elif isinstance(member, Strut):
                        vertex.set_strut(member)
                    else:
                        raise Exception('expected member to be instance of Strut or Tendon')
            print('*** Warning! xalglib_forces == None ***')

    @property
    def vertex_array(self):
        return array([vertex.coordinates for vertex in self.vertices])

    @property
    def strut_array(self):
        # s_array = []
        # for strut in self.struts:
        #     s_array.append([self.vertices.index(strut.vertices[0]), self.vertices.index(strut.vertices[1])])
        return [[self.vertices.index(strut.vertices[0]), self.vertices.index(strut.vertices[1])]
                for strut in self.struts]
        # return array(self.strut_vertices)

    @property
    def tendon_array(self):
        return [[self.vertices.index(tendon.vertices[0]), self.vertices.index(tendon.vertices[1])]
                for tendon in self.tendons]

    @property
    def xalglib_member_forces(self):
        return o1.solve_tensegrity_tensions(self.strut_array, self.tendon_array, self.vertex_array)

    @property
    def xalglib_vertex_forces(self):
        """ returns the net force on each vertex based on xalglib"""
        force_vector_list = []
        xalglib_forces = self.xalglib_member_forces
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
        return True not in [strut.spring_force_magnitude > err_tol for strut in self.struts]

    def print_spring_forces(self):
        type_width = len('vertex')
        coord_width = 20
        round_param = 2
        if self.equilibrium:
            print('*** This tensegrity is in equilibrium ***')
        else:
            print('*** This tensegrity is NOT in equilibrium ***')
        print('*      Coordinates          Tendon Force Vector')
        for vertex in self.vertices:
            print('Vertex',
                  f'{str(around(array(vertex.coordinates), round_param)): <{coord_width}}',
                  f'{str(around(array(vertex.tendon_spring_force_vector), round_param)): <{coord_width}}')
        print('*      V0 Coordinates       V1 Coordinates       Force      Unit Vector 0        Unit Vector 1')
        for tendon in self.tendons:
            # print(f'{0: <{type_width}}'.format(tendon.member_type),
            print(f'{tendon.member_type: <{type_width}}',
                  f'{str(around(array(tendon.vertices[0].coordinates), round_param)): <{coord_width}}',
                  f'{str(around(array(tendon.vertices[1].coordinates), round_param)): <{coord_width}}',
                  round(tendon.spring_force_magnitude, round_param), '     ',
                  f'{str(around(array(tendon.unit_force_vector(tendon.vertices[0])), round_param)): <{coord_width}}',
                  f'{str(around(array(tendon.unit_force_vector(tendon.vertices[1])), round_param)): <{coord_width}}'
                  )
        for strut in self.struts:
            print(f'{strut.member_type: <{type_width}}',
                  f'{str(around(array(strut.vertices[0].coordinates), round_param)): <{coord_width}}',
                  f'{str(around(array(strut.vertices[1].coordinates), round_param)): <{coord_width}}',
                  str(around(array(strut.spring_force_vector_sum), round_param)))


class Kite(Tensegrity):

    def __init__(self):
        coordinates = [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
        self.strut_vertices = [[0, 2], [1, 3]]
        self.tendon_vertices = [[0, 1], [1, 2], [2, 3], [3, 0]]
        Tensegrity.__init__(self, coordinates, self.strut_vertices, self.tendon_vertices)

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
            print('solver: spring_force_vector_sum',  self.struts[0].spring_force_vector_sum)
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

    def solver_strut(self):
        """ attempts to reach equilibrium by changing position of the top of the struts. Keeps all radial symmetries"""
        err_tol = 0.1
        max_steps = 1000
        # step_sizes = [5 * err_tol, 2 * err_tol, 0.9 * err_tol]
        step_size = 1.5 * err_tol  # 1.9 should work
        phi_direction = -1
        # i_step_size = 0
        step_count = 0
        while not self.equilibrium(err_tol=err_tol):
            # modify phi (spherical polar angle) of top vertex in all struts equally
            for strut in self.struts:
                # strut.modify_phi(phi_direction * step_sizes[i_step_size])
                initial_strut_force_mag = self.struts[0].spring_force_magnitude
                strut.modify_phi(phi_direction * step_size)
                if self.struts[0].spring_force_magnitude > initial_strut_force_mag:
                    # if we are going the wrong way then reverse direction
                    phi_direction = -phi_direction
                    # phi_better = False
                elif self.struts[0].spring_force_magnitude == initial_strut_force_mag:
                    raise Exception('no change in spring force found after change in phi')
                initial_strut_force_mag = self.struts[0].spring_force_magnitude
                strut.modify_phi(phi_direction * step_size)

                # else:
                # phi_better = True
            step_count += 1
            if step_count > max_steps:
                raise Exception('max step count exceeded')
            # modify theta (spherical azimuthal angle) of top vertex in all struts equally

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
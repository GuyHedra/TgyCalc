import math
import olof_1 as o1
from numpy import array

""" simpleTools.py is a simplified version of polyTools intended for use within TgyCalc"""


def add_vectors(v0, v1):
    return [c0 + c1 for c0, c1 in zip(v0, v1)]


def dot_product(v0, v1):
    """returns the dot product of two vectors"""
    return sum([v0_coord * v1_coord for v0_coord, v1_coord in zip(v0, v1)])


def normalize_vector(vector):
    magnitude = sum([element ** 2 for element in vector]) ** 0.5
    return [element / magnitude for element in vector]


class Vertex:

    def __init__(self, coordinates):
        self.coordinates = coordinates
        self.struts = []
        self.tendons = []

    def set_coordinates(self, coordinates):
        self.coordinates = coordinates

    def add_strut(self, strut):
        self.struts.append(strut)

    def add_tendon(self, tendon):
        self.tendons.append(tendon)

    @property
    def members(self):
        return self.struts + self.tendons

    def distance(self, other):
        return sum(
            [(p0_coord - p1_coord) ** 2 for p0_coord, p1_coord in zip(self.coordinates, other.coordinates)]) ** 0.5

    def __sub__(self, other):
        return [self_coord - other_coord for self_coord, other_coord in zip(self.coordinates, other.coordinates)]


class Member:
    """ base class for tendons and struts"""
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
    def force_magnitude(self):
        if self.current_length > self.nom_length:
            mag = self.spring_constant * (self.current_length - self.nom_length)
        else:
            mag = 0
        return mag

    def unit_force_vector(self, vertex):
        if (vertex is self.vertices[0] and isinstance(self, Strut)) or \
                (vertex is self.vertices[1] and isinstance(self, Tendon)):
            unit_vector = normalize_vector(self.vertices[0] - self.vertices[1])
        elif (vertex is self.vertices[1] and isinstance(self, Strut)) or \
                (vertex is self.vertices[0] and isinstance(self, Tendon)):
            unit_vector = normalize_vector(self.vertices[1] - self.vertices[0])
        else:
            raise Exception('Expected vertex to belong to member')
        return unit_vector


class Strut(Member):

    def __init__(self, vertices):
        Member.__init__(self, vertices)

    def unit_force_vector(self, vertex):
        if vertex is self.vertices[0]:
            unit_vector = normalize_vector(self.vertices[0] - self.vertices[1])
        elif vertex is self.vertices[1]:
            unit_vector = normalize_vector(self.vertices[1] - self.vertices[0])
        else:
            raise Exception('Expected vertex to belong to member')
        return unit_vector


class Tendon(Member):
    def __init__(self, vertices):
        Member.__init__(self, vertices)

    def unit_force_vector(self, vertex):
        if vertex is self.vertices[1]:
            unit_vector = normalize_vector(self.vertices[0] - self.vertices[1])
        elif vertex is self.vertices[0]:
            unit_vector = normalize_vector(self.vertices[1] - self.vertices[0])
        else:
            raise Exception('Expected vertex to belong to member')
        return unit_vector


class Tensegrity:

    def __init__(self, coordinates, strut_vertices, tendon_vertices):
        self.vertices = [Vertex(coords) for coords in coordinates]
        self.struts = [Strut([self.vertices[vtx_indices[0]], self.vertices[vtx_indices[1]]])
                       for vtx_indices in strut_vertices]
        self.tendons = [Tendon([self.vertices[vtx_indices[0]], self.vertices[vtx_indices[1]]])
                        for vtx_indices in tendon_vertices]
        self.populate_members()
        pass

    @property
    def members(self):
        return self.struts + self.tendons

    def populate_members(self):
        """ populate each vertex's list of tendons and struts (members) and store xalglib forces in each member"""
        xalglib_forces = self.xalglib_member_forces
        if len(xalglib_forces) != len(self.members):
            raise Exception('expected xalglib_forces to be same length as self.members')
        for member, force in zip(self.members, xalglib_forces):
            member.set_xalglib_force(force)
            for vertex in member.vertices:
                if isinstance(member, Tendon):
                    vertex.add_tendon(member)
                elif isinstance(member, Strut):
                    vertex.add_strut(member)
                else:
                    raise Exception('expected member to be instance of Strut or Tendon')

    @property
    def vertex_array(self):
        return array([vertex.coordinates for vertex in self.vertices])

    @property
    def strut_array(self):
        return array(self.strut_vertices)

    @property
    def tendon_array(self):
        return array(self.tendon_vertices)

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
                force_vector = add_vectors(force_vector, member.xalglib_force_vector(vertex))
                pass
            i_xalglib += 1
            force_vector_list.append(force_vector)
        return force_vector_list


class Kite(Tensegrity):

    def __init__(self):
        coordinates = [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
        self.strut_vertices = [[0, 2], [1, 3]]
        self.tendon_vertices = [[0, 1], [1, 2], [2, 3], [3, 0]]
        Tensegrity.__init__(self, coordinates, self.strut_vertices, self.tendon_vertices)


class PinnedKite(Tensegrity):
    """ struts are attached by a pivot at their midpoints and the angle between struts (theta) is changed to reach
    a force equilibrium """

    def __init__(self):
        coordinates = [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
        self.strut_vertices = [[0, 2], [1, 3]]
        self.tendon_vertices = [[0, 1], [1, 2], [2, 3], [3, 0]]
        Tensegrity.__init__(self, coordinates, self.strut_vertices, self.tendon_vertices)

    def set_strut_theta(self, theta):
        """ assumes that both struts have their midpoints at [0, 0, 0] and have a length of 2"""
        self.struts[0].vertices[0].set_coordinates([math.cos(theta), math.sin(theta), 0])
        self.struts[0].vertices[1].set_coordinates([-math.cos(theta), -math.sin(theta), 0])


class Prism(Tensegrity):

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
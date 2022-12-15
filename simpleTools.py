import math
import olof_1 as o1
from numpy import array

""" simpleTools.py is a simplified version of polyTools intended for use within TgyCalc"""


class Vertex:

    def __init__(self, coordinates):
        self.coordinates = coordinates
        self.struts = []
        self.tendons = []

    def add_strut(self, strut):
        self.struts.append(strut)

    def add_tendon(self, tendon):
        self.tendons.append(tendon)

    @property
    def members(self):
        return self.struts + self.tendons

    def __sub__(self, other):
        return [self_coord - other_coord for self_coord, other_coord in zip(self.coordinates, other.coordinates)]


class Strut:

    def __init__(self, vertices):
        self.vertices = vertices
        self.force = 0

    def force_vector(self, vertex):
        if vertex is self.vertices[0]:
            unit_vector = self.vertices[0] - self.vertices[1]
        elif vertex is self.vertices[1]:
            unit_vector = self.vertices[1] - self.vertices[0]
        else:
            raise Exception('Expected vertex to belong to member')
        return [self.force * coord for coord in unit_vector]


class Tendon:

    def __init__(self, vertices):
        self.vertices = vertices
        self.force = 0

    def force_vector(self, vertex):
        if vertex is self.vertices[0]:
            unit_vector = self.vertices[0] - self.vertices[1]
        elif vertex is self.vertices[1]:
            unit_vector = self.vertices[1] - self.vertices[0]
        return [self.force * coord for coord in unit_vector]


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
        """ populate each vertex's list of tendons and struts"""
        for member in self.members:
            for vertex in member.vertices:
                if member is Tendon:
                    vertex.add_tendon(member)
                elif member is Strut:
                    vertex.add_strut(member)

    @property
    def vertex_array(self):
        return array(self.coordinates)

    @property
    def strut_array(self):
        return array(self.strut_vertices)

    @property
    def tendon_array(self):
        return array(self.tendon_vertices)

    @property
    def member_forces(self):
        return o1.solve_tensegrity_tensions(self.strut_array, self.tendon_array, self.vertex_array)

    @property
    def vertex_forces(self):
        """ returns the net force on each vertex"""
        force_vector_list = []
        for vertex in self.vertices:
            force_vector = array([0, 0, 0])
            for member in vertex.members:
                member_vector = member.force_vector(vertex)
                force_vector = force_vector + array(member_vector)
                # force_vector = force_vector + array(member.force_vector(vertex))
            force_vector_list.append(force_vector)
        return force_vector_list


class Kite(Tensegrity):

    def __init__(self):
        self.coordinates = [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
        self.strut_vertices = [[0, 2], [1, 3]]
        self.tendon_vertices = [[0, 1], [1, 2], [2, 3], [3, 0]]
        Tensegrity.__init__(self, self.coordinates, self.strut_vertices, self.tendon_vertices)


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
""" forces.p implements a set of tools to calculate the forces acting on a tensegrity structure"""

import numpy as np


def other_vertex_index(member, vertex_index):
    if member[0] == vertex_index:
        return member[1]
    elif member[1] == vertex_index:
        return member[0]
    else:
        raise Exception('vertex_index not found in member')

def matrix_demo():
    pass


class Tensegrity:
    def __init__(self, vertices, struts, tendons):
        self.vertices = vertices
        self.struts = struts
        self.tendons = tendons
        self.tendon_force = 1
        self.strut_f_unit_vectors = {}  # key is vertex index: value is scalar strut force
        self.tendon_f_unit_vectors = {}  # key is vertex index: value is list of tendon forces
        """ Build dictionaries that let us find the struts and tendons that are connected to a particular vertex"""
        self.vertex_struts = {}  # key is vertex index: value is strut index
        self.vertex_tendons = {}  # key is vertex index: values is list of tendon indices
        for i, strut in enumerate(self.struts):
            if strut[0] in self.vertex_struts or strut[1] in self.vertex_struts:
                raise Exception('vertex has more than one associated strut')
            else:
                self.vertex_struts[strut[0]] = i
                self.vertex_struts[strut[1]] = i
        for i, tendon in enumerate(self.tendons):
            if tendon[0] in self.vertex_tendons:
                deb0 = self.vertex_tendons[tendon[0]]
                self.vertex_tendons[tendon[0]].extend([i])
                # self.vertex_tendons[tendon[0]] =
            else:
                self.vertex_tendons[tendon[0]] = [i]
            if tendon[1] in self.vertex_tendons:
                self.vertex_tendons[tendon[1]].extend([i])
            else:
                self.vertex_tendons[tendon[1]] = [i]

    def tendon_forces(self):
        """ We will use the convention that strut forces are positive and tendon forces are negative:
        For each vertex the sum of all forces is equal to 0 for a stable balanced tensegrity """
        """ Find the strut force vectors """
        for strut in self.struts:
            for i, vertex_index in enumerate(strut):
                # todo add check for multiple struts on a vertex
                if vertex_index in self.strut_f_unit_vectors:
                    raise Exception('multiple struts are associated with vertex', vertex_index)
                vector = self.vertices[vertex_index] - self.vertices[strut[(i + 1) % 2]]
                self.strut_f_unit_vectors[vertex_index] = vector/np.linalg.norm(vector)
        # print(self.strut_f_unit_vectors)
        """ Find the tendon force vectors """
        for tendon in self.tendons:
            for i, vertex_index in enumerate(tendon):
                vector = self.vertices[tendon[(i + 1) % 2]] - self.vertices[vertex_index]
                if vertex_index in self.tendon_f_unit_vectors:
                    self.tendon_f_unit_vectors[vertex_index].extend([vector/np.linalg.norm(vector)])
                else:
                    self.tendon_f_unit_vectors[vertex_index] = [vector/np.linalg.norm(vector)]
        print(self.tendon_f_unit_vectors)
        """ Build the force matrix to solve"""
        # we will start with vertex 0
        # forces = [self.tendon_f_unit_vectors[0] + self.strut_f_unit_vectors[0]]
        a = np.column_stack(self.tendon_f_unit_vectors[0])
        print('a', a)

        """ Solve the force matrix"""


if __name__ == '__main__':
    vertex_list = np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]])
    strut_list = [[0, 2], [1, 3]]
    tendon_list = [[0, 1], [1, 2], [2, 3], [3, 0]]
    matrix_demo()
    kite = Tensegrity(vertex_list, strut_list, tendon_list)
    kite.tendon_forces()

""" forces.p implements a set of tools to calculate the forces acting on a tensegrity structure"""

import numpy as np
import math


def other_vertex_index(member, vertex_index):
    if member[0] == vertex_index:
        return member[1]
    elif member[1] == vertex_index:
        return member[0]
    else:
        raise Exception('vertex_index not found in member')


def kite():
    vertex_list = np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]])
    strut_list = [[0, 2], [1, 3]]
    tendon_list = [[0, 1], [1, 2], [2, 3], [3, 0]]
    return vertex_list, strut_list, tendon_list


def prism(alpha=5*math.pi/6):
    vertices = np.array([[math.cos(0), math.sin(0), 0],
                      [math.cos(2*math.pi/3), math.sin(2*math.pi/3), 0],
                      [math.cos(4 * math.pi / 3), math.sin(4 * math.pi / 3), 0],
                      [math.cos(alpha), math.sin(alpha), 1],
                      [math.cos(alpha + 2 * math.pi / 3), math.sin(alpha + 2 * math.pi / 3), 1],
                      [math.cos(alpha + 4 * math.pi / 3), math.sin(alpha + 4 * math.pi / 3), 1]])
    #these numbers are indicies in the vertex_list array for the vertex_list on either end of the member
    compression_members = np.array([[0, 3],
                                 [1, 4],
                                 [2, 5]])
    tension_members = np.array([[0, 1],
                             [1, 2],
                             [2, 0],
                             [3, 4],
                             [4, 5],
                             [5, 3],
                             [1, 3],
                             [2, 4],
                             [0, 5]])
    return vertices, compression_members, tension_members


class Tensegrity:
    def __init__(self, vertices, struts, tendons, name):
        self.vertices = vertices
        self.struts = struts
        self.tendons = tendons
        self.name = name
        self.tendon_force = 1
        self.f_unit_vectors = {}  # key is vertex index: value is scalar strut force
        # self.f_unit_vectors = {}  # key is vertex index: value is list of tendon forces
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
            else:
                self.vertex_tendons[tendon[0]] = [i]
            if tendon[1] in self.vertex_tendons:
                self.vertex_tendons[tendon[1]].extend([i])
            else:
                self.vertex_tendons[tendon[1]] = [i]
        """ We will use the convention that strut forces are positive and tendon forces are negative:
        For each vertex the sum of all forces is equal to 0 for a stable balanced tensegrity """
        """ Find the strut force vectors """
        for strut in self.struts:
            for i, vertex_index in enumerate(strut):
                if vertex_index in self.f_unit_vectors:
                    raise Exception('multiple struts are associated with vertex', vertex_index)
                vector = self.vertices[vertex_index] - self.vertices[strut[(i + 1) % 2]]
                self.f_unit_vectors[vertex_index] = [vector / np.linalg.norm(vector)]
        # print(self.f_unit_vectors)
        """ Find the tendon force vectors """
        for tendon in self.tendons:
            for i, vertex_index in enumerate(tendon):
                vector = self.vertices[tendon[(i + 1) % 2]] - self.vertices[vertex_index]
                self.f_unit_vectors[vertex_index].extend([vector/np.linalg.norm(vector)])
        print(self.f_unit_vectors)

    def tendon_forces(self):
        """ Build the force matrix to solve"""
        # we will start with vertex 0
        # forces = [self.f_unit_vectors[0] + self.f_unit_vectors[0]]
        # a = np.column_stack([self.f_unit_vectors[0], [1, 0, 0]])
        a = np.column_stack([np.transpose(self.f_unit_vectors[0])])
        if self.name == 'kite':
            b = np.array([[0], [0], [1]])
            a[2] = [1, 0, 0]
        elif self.name == 'prism':
            b = np.array([[0], [0], [0], [1]])
            a = np.vstack((a, [1, 0, 0, 0]))
            # a.append = [1, 0, 0, 0]
        print('a', a)
        # print('a[2]', a[2])
        print('b', b)
        x = np.linalg.solve(a, b)
        # x = np.linalg.lstsq(a, b)
        print('x', x)
        """ Solve the force matrix"""


if __name__ == '__main__':
    # mode = 'kite'
    mode = 'prism'
    if mode == 'kite':
        thing = Tensegrity(*kite(), mode)
        thing.tendon_forces()
    elif mode == 'prism':
        thing = Tensegrity(*prism(), mode)
        thing.tendon_forces()

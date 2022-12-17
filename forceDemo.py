import math
from numpy import concatenate, shape, array, transpose, linalg, zeros, ones, full
# from xalglib import xalglib
import simpleTools as sT

""" demonstrate a method of adding the forces on a tensegrity structures vertices"""


if __name__ == '__main__':
    # print('*** kite forces ***')
    # kite = sT.Kite()
    # print('kite xalglib member forces', kite.xalglib_member_forces)
    # print('kite xalglib vertex forces', kite.xalglib_vertex_forces)
    # #
    # prism = sT.Prism()
    # print('*** prism forces')
    # print('prism xalglib member forces', prism.xalglib_member_forces)
    # print('prism xalglib vertex forces', prism.xalglib_vertex_forces)
    # #
    print('*** initial pinned kite forces ***')
    pinned_kite = sT.PinnedKite()
    pinned_kite.print_spring_forces()
    print('kite xalglib member forces', pinned_kite.xalglib_member_forces)
    print('kite xalglib vertex forces', pinned_kite.xalglib_vertex_forces)
    theta = math.pi/3
    print('\n*** twisted pinned kite forces ***')
    print('theta = ', theta)
    # pinned_kite.set_strut_theta(theta)
    pinned_kite.vertices[0].coordinates = [1.1, 0, 0]
    pinned_kite.print_spring_forces()
    print('kite xalglib member forces', pinned_kite.xalglib_member_forces)
    print('kite xalglib vertex forces', pinned_kite.xalglib_vertex_forces)

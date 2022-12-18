import math
from numpy import concatenate, shape, array, transpose, linalg, zeros, ones, full
# from xalglib import xalglib
import simpleTools as sT

""" demonstrate a method of adding the forces on a tensegrity structures vertices"""


if __name__ == '__main__':
    mode = 'standing prism solver'
    # print('*** kite forces ***')
    # kite = sT.Kite()
    # print('kite xalglib member forces', kite.xalglib_member_forces)
    # print('kite xalglib vertex forces', kite.xalglib_vertex_forces)
    #
    # prism = sT.PrismNeq3()
    # print('*** prism forces')
    # print('prism xalglib member forces', prism.xalglib_member_forces)
    # print('prism xalglib vertex forces', prism.xalglib_vertex_forces)
    #
    print('*** initial kite forces ***')
    # kite = sT.PinnedKite()
    kite = sT.Kite()
    kite.print_spring_forces()
    print('kite xalglib member forces', kite.xalglib_member_forces)
    # print('kite xalglib vertex forces', kite.xalglib_vertex_forces)
    theta = math.pi/3
    print('\n*** twisted pinned kite forces ***')
    print('theta = ', theta)
    # kite.set_strut_theta(theta)
    kite.vertices[0].coordinates = [2, 0, 0]
    kite.print_spring_forces()
    print('kite xalglib member forces', kite.xalglib_member_forces)
    print('kite xalglib vertex forces', kite.xalglib_vertex_forces)
    kite.solver_strut_0x()
    kite.print_spring_forces()
    # if mode == 'standing prism solver':
    #     prism = sT.Prism(n=3)
    #     prism.set_nom_tendon_lengths(0.8)
    #     prism.print_spring_forces()

import math
from numpy import concatenate, shape, array, transpose, linalg, zeros, ones, full
# from xalglib import xalglib
import TgyTools as sT

""" demonstrate a method of adding the forces on a tensegrity structures vertices"""


if __name__ == '__main__':
    mode = 'standing prism solver'
    # print('*** kite forces ***')
    # kite = TT.Kite()
    # print('kite xalglib member forces', kite.xalglib_member_forces)
    # print('kite xalglib terminal_vertex forces', kite.xalglib_vertex_forces)
    #
    # prism = TT.PrismNeq3()
    # print('*** prism forces')
    # print('prism xalglib member forces', prism.xalglib_member_forces)
    # print('prism xalglib terminal_vertex forces', prism.xalglib_vertex_forces)
    #
    # print('*** initial kite forces ***')
    # kite = TT.PinnedKite()
    print('* Initializing Kite *')
    kite = sT.Kite()
    kite.print_spring_forces()
    # print('kite xalglib member forces', kite.xalglib_member_forces)
    # # print('kite xalglib terminal_vertex forces', kite.xalglib_vertex_forces)
    # theta = math.pi/3
    # print('\n*** twisted pinned kite forces ***')
    # print('theta = ', theta)
    # kite.set_strut_theta(theta)
    print('*** modifying kite vertices such that it will be out of equilibrium')
    kite.vertices[0].coordinates = [2, 0, 0]
    kite.print_spring_forces()
    print('*** updating xalglib forces')
    kite.set_xalglib_forces()
    print('kite xalglib member forces', kite.xalglib_member_forces)
    print('kite xalglib terminal_vertex forces', kite.xalglib_vertex_forces)
    print('* solving kite *')
    kite.solver_strut_0x()
    kite.print_spring_forces()
    print('*** updating xalglib forces')
    kite.set_xalglib_forces()
    print('kite xalglib member forces', kite.xalglib_member_forces)
    print('kite xalglib terminal_vertex forces', kite.xalglib_vertex_forces)
    # if mode == 'standing prism solver':
    #     prism = TT.Prism(n=3)
    #     prism.set_nom_tendon_lengths(0.8)
    #     prism.print_spring_forces()

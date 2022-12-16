import math
from numpy import concatenate, shape, array, transpose, linalg, zeros, ones, full
# from xalglib import xalglib
import simpleTools as sT

""" demonstrate a method of adding the forces on a tensegrity structures vertices"""


if __name__ == '__main__':
    print('*** kite forces ***')
    kite = sT.Kite()
    print('kite xalglib member forces', kite.xalglib_member_forces)
    print('kite xalglib vertex forces', kite.xalglib_vertex_forces)
    #
    prism = sT.Prism()
    print('*** prism forces')
    print('prism xalglib member forces', prism.xalglib_member_forces)
    print('prism xalglib vertex forces', prism.xalglib_vertex_forces)
    #
    pinned_kite = sT.PinnedKite()
    print('*** initial pinned kite forces ***')
    print('kite xalglib member forces', pinned_kite.xalglib_member_forces)
    print('kite xalglib vertex forces', pinned_kite.xalglib_vertex_forces)
    pinned_kite.set_strut_theta(math.pi/10)
    print('*** twisted pinned kite forces ***')
    print('kite xalglib member forces', pinned_kite.xalglib_member_forces)
    print('kite xalglib vertex forces', pinned_kite.xalglib_vertex_forces)

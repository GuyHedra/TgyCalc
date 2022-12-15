import math
from numpy import concatenate, shape, array, transpose, linalg, zeros, ones, full
# from xalglib import xalglib
import simpleTools as sT

""" demonstrate a method of adding the forces on a tensegrity structures vertices"""


if __name__ == '__main__':
    print('*** kite forces ***')
    kite = sT.Kite()
    print('kite member forces', kite.member_forces)
    print('kite vertex forces', kite.vertex_forces)
    #
    prism = sT.Prism()
    print('*** prism forces')
    print('prism member forces', prism.member_forces)
    print('prism vertex forces', prism.vertex_forces)

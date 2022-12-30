from scipy.spatial.transform import Rotation as R
import math
from numpy import array, set_printoptions
""" demonstrate a vector rotation using scipy"""


def rotation(angle):
    r1 = R.from_euler('z', angle)
    r2 = R.from_rotvec(angle * array([0, 0, 1]))
    v0 = [1, 0, 0]
    v1 = v0
    v2 = v1
    print('v0', v0)
    print('v1', r1.apply(v1))
    print('v2', r2.apply(v2))


def numpy_array_demo():
    set_printoptions(formatter={'float': '{: 10.3f}'.format})
    a = array([1, 2, 3])
    b = array([4, 5, 6])
    print('a + b', a + b)
    print('array 10.3f 123456789012345678901234567890')
    print('array 10.3f', array([-123456789.01, -123456789.0]))
    print('array 10.3f', array([-12345678.01, -12345678.0]))
    print('array 10.3f 123456789012345678901234567890')
    print('array 10.3f', array([-1234567.01, -1234567.0]))
    print('array 10.3f 123456789012345678901234567890')
    print('array 10.3f', array([-123456.01, -123456.0]))
    print('array 10.3f 123456789012345678901234567890')
    print('array 10.3f', array([-12345.01, -12345.0]))
    print('array 10.3f 123456789012345678901234567890')
    print('array 10.3f', array([-1234.01, -1234.0]))



class Vertex:

    def __init__(self, coordinates):
        self.coordinates = array(coordinates)


def vertex_demo():
    v0 = Vertex([1, 2, 3])
    print('v0', v0.coordinates)
    v_offset = array([4, 5, 6])
    print('v0 + v_offset', v0.coordinates + v_offset)
    print('v0 * v_offset', v0.coordinates * v_offset)
    print('v0 * 7', v0.coordinates * 7)


if __name__ == '__main__':
    angle = math.pi/4
    # rotation(angle)
    numpy_array_demo()
    # vertex_demo()

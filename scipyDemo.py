from scipy.spatial.transform import Rotation as R
import math
from numpy import array
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

if __name__ == '__main__':
    angle = math.pi/4
    rotation(angle)

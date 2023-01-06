import unittest
import math
import TgyPar as typ
import numpy as np
from numpy.testing import assert_almost_equal


class TestTrilateration(unittest.TestCase):
    def test_known_values(self):
        # [p1, p2, p3, r1, r2, r3, x1, x2]
        known_values = [
            [[(8 / 9) ** 0.5, 0, -1 / 3],
             [-(2 / 9) ** 0.5, (2 / 3) ** 0.5, -1 / 3],
             [-(2 / 9) ** 0.5, -(2 / 3) ** 0.5, -1 / 3],
             (8 / 3) ** 0.5, (8 / 3) ** 0.5, (8 / 3) ** 0.5,
             [0, 0, 1], [0, 0, -5 / 3]],
            [[1, 0, -1 / 2 ** 0.5],
             [-1, 0, -1 / 2 ** 0.5],
             [0, 1, 1 / 2 ** 0.5],
             2, 2, 2,
             [0, 5 / 3, -1.1785113019775793],
             [0, -1, 1 / 2 ** 0.5]]
        ]
        for p1, p2, p3, r1, r2, r3, x1, x2 in known_values:
            k1, k2 = typ.trilateration(p1, p2, p3, r1, r2, r3)
            assert_almost_equal(np.array(k1), np.array(x1))
            assert_almost_equal(np.array(k2), np.array(x2))


class TestDistance(unittest.TestCase):
    def test_distance(self):
        known_values = [[[1, 0, 0], [0, 0, 0], 1],
                        [[1, 1, 0], [-1, -1, 0], 2 ** 0.5 * 2],
                        [[2, 2, 0], [0, 0, 0], 2 ** 0.5 * 2]]
        for p0, p1, x in known_values:
            # v0 = typ.Vertex(c0, level=0)
            # v1 = typ.Vertex(c1, level=0)
            # print('c0, c1, x, d', c0, c1, x, v0.distance(v1))
            self.assertAlmostEqual(typ.distance(p0, p1), x)

class TestVertex(unittest.TestCase):
    def test_distance(self):
        known_values = [[[1, 0, 0], [0, 0, 0], 1],
                        [[1, 1, 0], [-1, -1, 0], 2 ** 0.5 * 2],
                        [[2, 2, 0], [0, 0, 0], 2 ** 0.5 * 2]]
        for c0, c1, x in known_values:
            v0 = typ.Vertex(c0, level=0)
            v1 = typ.Vertex(c1, level=0)
            # print('c0, c1, x, d', c0, c1, x, v0.distance(v1))
            self.assertAlmostEqual(v0.distance(v1), x)


if __name__ == '__main__':
    unittest.main()

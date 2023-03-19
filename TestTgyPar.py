import unittest
import math
import TgyPar as TP
import numpy as np
from numpy.testing import assert_almost_equal


class PrismRef:
    def __init__(self, name):
        if name == 'n=3, levels=1':
            """ From Geodesic Math and how to use it"""
            self.n = 3
            self.heights = [1]
            self.radii = [0.5, 0.5]
            self.strut_length = 1.3903
            self.bot_top = 1.03
            self.waist = 0.866
        else:
            raise Exception('name', name, 'is not supported')


class BoxRef:
    """ A reference tensegrity with all angles = pi/4 for testing only. Not stable."""
    top_coords = [[1, 1, 1], [-1, 1, 1], [-1, -1, 1], [1, -1, 1]]
    bot_coords = [[1, 1, -1], [-1, 1, -1], [-1, -1, -1], [1, -1, -1]]
    top_vertices = [TP.Vertex(np.array(coords), level=0) for coords in top_coords]
    bot_vertices = [TP.Vertex(np.array(coords), level=0) for coords in bot_coords]


class TestTowerPrism(unittest.TestCase):
    def test_known_values(self):
        """ specify strut length"""
        ref_tower = PrismRef('n=3, levels=1')
        tower = TP.PrismTower(n=ref_tower.n, levels=1, radii=ref_tower.radii, strut_lengths=[ref_tower.strut_length])
        assert_almost_equal(tower.struts[0].curr_length, ref_tower.strut_length, 4)
        tower.print_build()
        """ specify height """
        ref_tower = PrismRef('n=3, levels=1')
        tower = TP.PrismTower(n=ref_tower.n, levels=1, radii=ref_tower.radii, heights=ref_tower.heights)
        assert_almost_equal(tower.struts[0].curr_length, ref_tower.strut_length, 4)


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
            k1, k2 = TP.trilateration(p1, p2, p3, r1, r2, r3)
            assert_almost_equal(np.array(k1), np.array(x1))
            assert_almost_equal(np.array(k2), np.array(x2))


class TestDistance(unittest.TestCase):
    def test_distance(self):
        known_values = [[[1, 0, 0], [0, 0, 0], 1],
                        [[1, 1, 0], [-1, -1, 0], 2 ** 0.5 * 2],
                        [[2, 2, 0], [0, 0, 0], 2 ** 0.5 * 2]]
        for p0, p1, x in known_values:
            # v0 = TP.Vertex(c0, level=0)
            # v1 = TP.Vertex(c1, level=0)
            # print('c0, c1, x, d', c0, c1, x, v0.distance(v1))
            self.assertAlmostEqual(TP.distance(p0, p1), x)


class TestVertex(unittest.TestCase):
    def test_distance(self):
        known_values = [[[1, 0, 0], [0, 0, 0], 1],
                        [[1, 1, 0], [-1, -1, 0], 2 ** 0.5 * 2],
                        [[2, 2, 0], [0, 0, 0], 2 ** 0.5 * 2]]
        for c0, c1, x in known_values:
            v0 = TP.Vertex(c0, level=0)
            v1 = TP.Vertex(c1, level=0)
            # print('c0, c1, x, d', c0, c1, x, v0.distance(v1))
            self.assertAlmostEqual(v0.distance(v1), x)


class TestAngle(unittest.TestCase):
    def test_known_values(self):
        p0 = [0, 0, 1]
        p1 = [1, 0, 0]
        value = TP.angle(p0, p1)
        self.assertAlmostEqual(value, math.pi / 2)


if __name__ == '__main__':
    unittest.main()

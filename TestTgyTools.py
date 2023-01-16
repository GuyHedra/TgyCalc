import unittest
import TgyTools as TT
import math
import numpy as np
from numpy.testing import assert_almost_equal


class TgyTools(unittest.TestCase):
    # def test_free_vector(self):
    #     f_vector = TT.Vertex([1, 0, 0])
    #     self.assertEqual(f_vector.magnitude, 1)
    #     r_vector = f_vector.rotate(axis=[0, 0, 1], angle=math.pi/4)
    #     self.assertAlmostEqual(r_vector[0], 2 ** 0.5 / 2)

    def test_vertex(self):
        v1 = TT.Vertex([1, 0, 0])
        r1, theta1, z1 = v1.cyl_coordinates_rad
        self.assertEqual(r1, 1)
        self.assertEqual(theta1, 0)
        self.assertEqual(z1, 0)
        center = TT.Vertex([0, 0, 0])
        v1.rotate(center=center, axis=[0, 0, 1], angle=math.pi/4)
        # print('v1.coordinates', v1.coordinates)
        self.assertAlmostEqual(v1.coordinates[0], 2 ** 0.5 / 2)
        self.assertAlmostEqual(v1.coordinates[1], 2 ** 0.5 / 2)
        self.assertEqual(v1.coordinates[2], 0)
        r1_r, theta1_r, z1_r = v1.cyl_coordinates_rad
        self.assertEqual(r1_r, 1)
        self.assertAlmostEqual(theta1_r, math.pi/4)
        self.assertEqual(z1_r, 0)
        v2 = TT.Vertex([1, 0, 0])
        center = TT.Vertex([2, 0, 0])
        v2.rotate(center=center, axis=[0, 0, 1], angle=math.pi/4)
        # print('v2.coordinates', v2.coordinates)
        self.assertAlmostEqual(v2.coordinates[0], 2 - 2 ** 0.5 / 2)
        self.assertAlmostEqual(v2.coordinates[1], -2 ** 0.5 / 2)
        self.assertEqual(v2.coordinates[2], 0)
        v3 = TT.Vertex([-1, 1, 1])
        r3, theta3, z3 = v3.cyl_coordinates_rad
        self.assertAlmostEqual(r3, 2 ** 0.5)
        self.assertAlmostEqual(theta3, 3 * math.pi / 4)
        self.assertEqual(z3, 1)
        v4 = TT.Vertex([-1, -1, -1])
        r4, theta4, z4 = v4.cyl_coordinates_rad
        self.assertAlmostEqual(r4, 2 ** 0.5)
        self.assertAlmostEqual(theta4, -3 * math.pi / 4)
        self.assertEqual(z4, -1)
        v5 = TT.Vertex([1, -1, -2])
        r5, theta5, z5 = v5.cyl_coordinates_rad
        self.assertAlmostEqual(r5, 2 ** 0.5)
        self.assertAlmostEqual(theta5, -math.pi / 4)
        self.assertEqual(z5, -2)

    def test_vector_projection(self):
        # known values [unit_vector, proj_vector, x_vector]
        known_values = [[[1, 0, 0], [1, 1, 0], [1, 0, 0]],
                        [[2, 0, 0], [1, 1, 0], [1, 0, 0]],
                        [[1, 1, 1], [1, 0, 0], [1/3, 1/3, 1/3]],
                        ]
        for unit_vector, proj_vector, x_vector in known_values:
            assert_almost_equal(np.array(TT.vector_projection(proj_vector, unit_vector)), np.array(x_vector))

    def test_tendon_lateral_f_vec(self):
        coordinates = [[0, 1, 0], [0, -1, 0], [1, 1, 0], [1, 0, 0]]
        strut_vertices = [[0, 1]]
        tendon_vertices = [[0, 2], [1, 3]]
        nom_tendon_lengths = 0.5
        tensegrity = TT.TArbitrary(coordinates, strut_vertices, tendon_vertices, nom_tendon_lengths)
        vertex0 = tensegrity.vertices[0]
        vertex1 = tensegrity.vertices[1]
        tensegrity.tendons[0].set_force(1)
        tensegrity.tendons[1].set_force(1)
        assert_almost_equal(np.array(tensegrity.tendons[0].lateral_f_vec(vertex0)), np.array([1, 0, 0]))
        assert_almost_equal(np.array(tensegrity.tendons[1].lateral_f_vec(vertex1)), np.array([2 ** 0.5 / 2, 0, 0]))

    # def test_balance_forces(self):
    #     tgy = TT.Prism(n=3)
    #     tgy.balance_forces(verbose=2)

    def test_tensegrity_plot(self):
        tgy = TT.Prism(n=3)
        tgy.set_nom_tendon_lengths(0.1)
        tgy.balance_forces(verbose=2)
        tgy.plot(lateral_f=True, vertex_f=True)


if __name__ == '__main__':
    unittest.main()

import unittest
import simpleTools as sT
import math


class SimpleTools(unittest.TestCase):
    # def test_free_vector(self):
    #     f_vector = sT.Vertex([1, 0, 0])
    #     self.assertEqual(f_vector.magnitude, 1)
    #     r_vector = f_vector.rotate(axis=[0, 0, 1], angle=math.pi/4)
    #     self.assertAlmostEqual(r_vector[0], 2 ** 0.5 / 2)

    def test_vertex(self):
        v1 = sT.Vertex([1, 0, 0])
        v1.rotate(center=[0, 0, 0], axis=[0, 0, 1], angle=math.pi/4)
        print('v1.coordinates', v1.coordinates)
        self.assertAlmostEqual(v1.coordinates[0], 2 ** 0.5 / 2)
        self.assertAlmostEqual(v1.coordinates[1], 2 ** 0.5 / 2)
        self.assertEqual(v1.coordinates[2], 0)
        v2 = sT.Vertex([1, 0, 0])
        v2.rotate(center=[2, 0, 0], axis=[0, 0, 1], angle=math.pi/4)
        print('v2.coordinates', v2.coordinates)
        self.assertAlmostEqual(v2.coordinates[0], 2 - 2 ** 0.5 / 2)
        self.assertAlmostEqual(v2.coordinates[1], -2 ** 0.5 / 2)
        self.assertEqual(v2.coordinates[2], 0)



if __name__ == '__main__':
    unittest.main()

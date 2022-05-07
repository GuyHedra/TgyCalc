import unittest
import math
import pickle
import polyTools as pt
import things
# import testui


class test_plane_projection(unittest.TestCase):
    def testplaneprojectionknownvalues(self):
        # [point, normal, vector_v0, vector_v1, projection]
        cases = [[(0, 0, 0), (-1, 0, 1), (0, 0, 0), (10, 0, 0), (5, 0, 5)],
                 [(0, 0, 3), (0, 0, 1), (0, 0, 0), (10, 10, 10), (10, 10, 0)],
                 [(0, 0, 0), (1, 1, 1), (0, 0, 0), (0, 0, 10),
                  (-3.333333333333334, -3.333333333333334, 6.666666666666666)]
                 ]
        for case in cases:
            point_coords, normal_coords, vector_v0_coords, vector_v1_coords, projection_coords = case
            point = pt.Point(point_coords)
            normal = pt.Point(normal_coords)
            plane = pt.Plane(point, normal)
            vector = pt.Vector(pt.Point(vector_v0_coords), pt.Point(vector_v1_coords))
            projection = plane.projection(vector)
            print(projection, projection.magnitude)
            self.assertEqual(projection_coords, projection)



# class Test_Tendon(unittest.TestCase):
#     def testtendonfunctionality(self):
#         self.tensegrity_type = 'prism'
#         ph = pt.Tensegrity(self.tensegrity_type)
#         # strut = blend.add_strut_mesh(ph, strut_type=self.strut_type,
#         #                              strut_radii=self.strut_radii, strut_bend_radius=self.strut_bend_radius,
#         #                              strut_ring_count=self.strut_ring_count,
#         #                              strut_segment_count=self.strut_segment_count,
#         #                              strut_wall_thickness=self.strut_wall_thickness,
#         #                              vertex_group_name='strut')
#         # struts = blend.draw_struts(strut, ph)
#         # if self.strut_type in ['cylinder', 'tube']:
#         # #     ph.tendon_end_points(self.strut_radii)
#         # tendon = blend.add_tendon_mesh(tendon_radius=self.tendon_radius,
#         #                                tendon_segment_count=self.tendon_segment_count)
#         # ph.set_tendon_end_points((2,2,2))
#         # drawn_length = ph.tendons[0].drawn_length
#         main = ph.tendons[0].main
#         for tendon in ph.tendons:
#             for vertex, end_point_vector in zip(tendon.vertices, tendon.end_point_vectors):
#                 self.assertTrue(vertex in vertex.strut.vertices)
#                 self.assertTrue(vertex is end_point_vector.points[0])
#
#
# class TestVertex(unittest.TestCase):
#     def testvectorfunctionality(self):
#         test_vertex0 = pt.Vertex(pt.Point((1, 2, 3)))
#         test_vertex1 = pt.Vertex(pt.Point((4, 5, 6)))
#         my_coords = (test_vertex1 - test_vertex0).coords
#         # print('my_coords', my_coords)
#
#
# class TestPoint(unittest.TestCase):
#     def testpointrotate(self):
#         test_point = pt.Point((1, 0, 0))
#         rotated_point = test_point.rotate((0, 0, 0), pt.Point((0, 0, 1)), math.pi)
#         self.assertTrue(rotated_point == (-1, 0, 0))

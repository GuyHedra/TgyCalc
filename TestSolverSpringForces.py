import unittest
import simpleTools as sT
import math


class Solver(unittest.TestCase):
    # def test_x_displace_kite(self):
    #     print('>>>Initializing kite')
    #     kite = sT.Kite()
    #     kite.print_spring_forces()
    #     print('>>>Pushing strut 0 by 1 in x')
    #     displacement = [1, 0, 0]
    #     kite.struts[0].vertices[0].set_coordinates(sT.vector_add(kite.struts[0].vertices[0].coordinates, displacement))
    #     kite.struts[0].vertices[1].set_coordinates(sT.vector_add(kite.struts[0].vertices[1].coordinates, displacement))
    #     kite.print_spring_forces()
    #     kite.solver_spring_forces()
    #     kite.print_spring_forces()
    #     self.assertTrue(kite.equilibrium(err_tol=0.001))

    def test_theta_displace_kite(self):
    #     err_tol = 0.01
    #     verbose = False
    #     print('>>>Initializing kite')
    #     kite = sT.Kite()
    #     # anchor one vertex of each strut to make results easier to interpret
    #     kite.vertices[2].set_anchor(True)
    #     kite.vertices[3].set_anchor(True)
    #     kite.print_spring_forces(err_tol=err_tol, vertices=False)
    #     print('>>>Rotate strut 0 by pi/4')
    #     angle = math.pi / 4
    #     kite.struts[0].rotate(center_vertex=kite.struts[0].vertices[1], axis=[0, 0, 1], angle=angle)
    #     kite.print_spring_forces(err_tol=err_tol, vertices=False)
    #     kite.solver_spring_forces(err_tol=err_tol, verbose=verbose)
    #     kite.print_spring_forces(err_tol=err_tol, vertices=False)
    #     self.assertTrue(kite.equilibrium(err_tol=err_tol))
    #     # do it again without the anchors
    #     # verbose = True
    #     verbose = False
    #     print('>>>Initializing kite')
    #     kite = sT.Kite()
    #     # anchor one vertex of each strut to make results easier to interpret
    #     kite.print_spring_forces(err_tol=err_tol, vertices=False)
    #     print('>>>Rotate strut 0 by pi/4')
    #     angle = math.pi / 4
    #     kite.struts[0].rotate(center_vertex=kite.struts[0].vertices[1], axis=[0, 0, 1], angle=angle)
    #     kite.print_spring_forces(err_tol=err_tol, vertices=False)
    #     err_tol = 0.2
    #     kite.solver_spring_forces(err_tol=err_tol, max_step_count=1000, initial_step=err_tol/2, verbose=verbose)
    #     kite.print_spring_forces(err_tol=err_tol, vertices=False)
    #     self.assertTrue(kite.equilibrium(err_tol=err_tol))
        # do it again without the anchors with parameters that will cause an overshoot
        err_tol = 0.1
        verbose = True
        # verbose = False
        print('>>>Initializing kite')
        kite = sT.Kite()
        # anchor one vertex of each strut to make results easier to interpret
        kite.print_spring_forces(err_tol=err_tol, vertices=False)
        print('>>>Rotate strut 0 by pi/4')
        angle = math.pi / 4
        kite.struts[0].rotate(center_vertex=kite.struts[0].vertices[1], axis=[0, 0, 1], angle=angle)
        kite.print_spring_forces(err_tol=err_tol, vertices=False)
        kite.solver_spring_forces(err_tol=err_tol, max_step_count=1000, initial_step=err_tol/10, verbose=verbose)
        kite.print_spring_forces(err_tol=err_tol, vertices=False)
        self.assertTrue(kite.equilibrium(err_tol=err_tol))


if __name__ == '__main__':
    unittest.main()

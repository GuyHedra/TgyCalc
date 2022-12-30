import unittest
import TgyTools as TT
import olof_1 as o1
from numpy import array
import math


class Solver(unittest.TestCase):
    def test_x_displace_kite(self):
        err_tol = 0.1
        verbosity = 0
        verbosity = 2
        if verbosity > 0:
            print('>>>Initializing kite')
        kite = TT.Kite()
        kite.print_spring_forces(err_tol=err_tol, verbosity=verbosity)
        if verbosity > 0:
            print('>>>Pushing strut 0 by 1 in x')
        displacement = [1, 0, 0]
        kite.struts[0].vertices[0].set_coordinates(TT.vector_add(kite.struts[0].vertices[0].coordinates, displacement))
        kite.struts[0].vertices[1].set_coordinates(TT.vector_add(kite.struts[0].vertices[1].coordinates, displacement))
        kite.print_spring_forces(err_tol=err_tol, verbosity=verbosity)
        kite.solver()
        kite.print_spring_forces(err_tol=err_tol, verbosity=verbosity)
        if verbosity > 0:
            print('algblib ', o1.solve_tensegrity_tensions(array(kite.strut_array),
                                                           array(kite.tendon_array),
                                                           array(kite.vertex_array)))
        self.assertTrue(kite.equilibrium(err_tol=err_tol))

    def test_solver_spring_forces(self):
        """ test each of the 7 operation cases that solver_spring_forces uses """
        run_cases = [1, 2, 3, 4, 5, 6]
        run_cases = [4]
        # Both vertex forces lateral and more parallel than orthogonal:
        # Translate strut along sum of vertex forces
        # first make vertex forces completely parallel
        err_tol = 0.1
        verbosity = 0
        # verbosity = 2
        if 1 in run_cases:
        # case 1
            if verbosity > 0:
                print('****** test case 1 ******')
            coordinates = [[0, 1, 0], [0, -1, 0], [1, 1, 0], [1, -1, 0]]
            strut_vertices = [[0, 1]]
            tendon_vertices = [[0, 2], [1, 3]]
            nom_tendon_lengths = 0.5
            tensegrity = TT.TArbitrary(coordinates, strut_vertices, tendon_vertices, nom_tendon_lengths)
            tensegrity.print_spring_forces(err_tol=err_tol, verbosity=verbosity)
            tensegrity.solver(err_tol=err_tol, max_step_count=1000, initial_step=err_tol / 10, verbose=False)
            tensegrity.print_spring_forces(err_tol=err_tol, verbosity=verbosity)
            self.assertTrue(tensegrity.equilibrium(err_tol=err_tol))
        # move the anchor vertices apart in the y direction
        if 2 in run_cases:
            verbosity = 0
            # verbosity = 2
            if verbosity > 0:
                print('****** test case 2 ******')
            coordinates = [[0, 1, 0], [0, -1, 0], [1, 1.9, 0], [1, -1.9, 0]]
            strut_vertices = [[0, 1]]
            tendon_vertices = [[0, 2], [1, 3]]
            nom_tendon_lengths = 0.5
            tensegrity = TT.TArbitrary(coordinates, strut_vertices, tendon_vertices, nom_tendon_lengths)
            tensegrity.print_spring_forces(err_tol=err_tol, verbosity=verbosity)
            tensegrity.solver(err_tol=err_tol, max_step_count=1000, initial_step=err_tol / 10, verbose=False)
            tensegrity.print_spring_forces(err_tol=err_tol, verbosity=verbosity)
            self.assertTrue(tensegrity.equilibrium(err_tol=err_tol))
        # put one anchor .5 directly above top of strut and other anchor 1.0 to right of bot of strut
        if 3 in run_cases:
            verbosity = 0
            # verbosity = 2
            if verbosity > 0:
                print('****** test case 3 ******')
            coordinates = [[0, 1, 0], [0, -1, 0], [0, 1.5, 0], [1, -1, 0]]
            strut_vertices = [[0, 1]]
            tendon_vertices = [[0, 2], [1, 3]]
            nom_tendon_lengths = 0.4
            tensegrity = TT.TArbitrary(coordinates, strut_vertices, tendon_vertices, nom_tendon_lengths)
            tensegrity.print_spring_forces(err_tol=err_tol, verbosity=verbosity)
            tensegrity.solver(err_tol=err_tol, max_step_count=1000, initial_step=err_tol / 10, verbose=False)
            tensegrity.print_spring_forces(err_tol=err_tol, verbosity=verbosity)
            self.assertTrue(tensegrity.equilibrium(err_tol=err_tol))
        if 4 in run_cases:
            verbosity = 0
            verbosity = 2
            if verbosity > 0:
                print('****** test case 4 Prism with n = 3 ******')
            # tensegrity = TT.Prism(n=4)
            tensegrity = TT.Prism(n=3)
            err_tol = 1
            # if verbosity > 0:
            #     print('algblib ', o1.solve_tensegrity_tensions(array(tensegrity.strut_array),
            #                                                    array(tensegrity.tendon_array),
            #                                                    array(tensegrity.vertex_array)))
            # set waist tendon spring constants to .1 of vertical tendons spring constants
            for tendon in tensegrity.tendons[0:6]:
                tendon.set_spring_constant(0.1)
            tensegrity.set_nom_tendon_lengths(0.5)
            # tensegrity.print_spring_forces(err_tol=err_tol, verbosity=verbosity)
            tensegrity.print_cylindrical(err_tol=err_tol)
            if verbosity > 0:
                print('*** solving test case 4 ***')
            # tensegrity.solver(err_tol=err_tol, max_step_count=1000, initial_step=err_tol / 10, verbose=False)
            # debug

            tensegrity.solver(err_tol=err_tol, max_step_count=100, initial_step=err_tol / 10, verbose=True)
            tensegrity.print_cylindrical(err_tol=err_tol)
            verbosity = 0
            if verbosity > 0:
                print('algblib ', o1.solve_tensegrity_tensions(array(tensegrity.strut_array),
                                                               array(tensegrity.tendon_array),
                                                               array(tensegrity.vertex_array)))
            self.assertTrue(tensegrity.equilibrium(err_tol=err_tol))
        if 5 in run_cases:
            verbosity = 0
            if verbosity > 0:
                print('****** test case 5 ******')
            err_tol = 0.1
            coordinates = [[0, 1, 0], [0, -1, 0], [1, 0, 0], [-1, 0, 0]]
            strut_vertices = [[0, 1]]
            tendon_vertices = [[0, 2], [0, 3], [1, 2], [1, 3]]
            nom_tendon_lengths = 0.4
            tensegrity = TT.TArbitrary(coordinates, strut_vertices, tendon_vertices, nom_tendon_lengths)
            tensegrity.vertices[0].coordinates = [0.5, 1, 1]
            tensegrity.vertices[1].coordinates = [0.5, -1, -1]
            tensegrity.print_spring_forces(err_tol=err_tol, verbosity=verbosity)
            if verbosity > 0:
                print('****** solving test case 5 ******')
            tensegrity.solver(err_tol=err_tol, max_step_count=1000, initial_step=err_tol / 10, verbose=False)
            tensegrity.print_spring_forces(err_tol=err_tol, verbosity=verbosity)
            self.assertTrue(tensegrity.equilibrium(err_tol=err_tol))
        if 6 in run_cases:
            verbosity = 2
            if verbosity > 0:
                print('****** test case 6 ******')
            err_tol = 0.1
            coordinates = [[0, 1, 0], [0, -1, 0], [1, 0, 0], [-1, 0, 0]]
            strut_vertices = [[0, 1]]
            tendon_vertices = [[0, 2], [0, 3], [1, 2], [1, 3]]
            nom_tendon_lengths = 0.4
            tensegrity = TT.TArbitrary(coordinates, strut_vertices, tendon_vertices, nom_tendon_lengths)
            tensegrity.vertices[0].coordinates = [0.5, 1, 1]
            tensegrity.vertices[1].coordinates = [0.5, -1, -1]
            tensegrity.print_spring_forces(err_tol=err_tol, verbosity=verbosity)
            if verbosity > 0:
                print('****** solving test case 6 ******')
            tensegrity.solver(err_tol=err_tol, max_step_count=1000, initial_step=err_tol / 10, verbose=False)
            tensegrity.print_spring_forces(err_tol=err_tol, verbosity=verbosity)
            self.assertTrue(tensegrity.equilibrium(err_tol=err_tol))

if __name__ == '__main__':
    unittest.main()

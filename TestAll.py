import unittest
import TestTgyTools
import TestSolverSpringForces
import time

class TestTgyCalc(unittest.TestCase):
    def test_all(self):
        testsuite = unittest.defaultTestLoader.loadTestsFromModule(TestTgyTools)
        testsuite.addTests(unittest.defaultTestLoader.loadTestsFromModule(TestSolverSpringForces))
        runner = unittest.TextTestRunner(verbosity=2)
        runner.run(testsuite)

# if __name__ == '__main__':
#     unittest.main()
#
# import unittest
# import time
# # import profile
# import testpolyhedra, testthings, testph_types
#
# # todo: create a test script to run in Blender environment, preferably from the command line
# class TestPolyhedra(unittest.TestCase):
#     """ Run all of the tests in the Polyhedra project"""
#     def testpolyhedra(self):
#         testsuite = unittest.defaultTestLoader.loadTestsFromModule(testthings)
#         testsuite.addTests(unittest.defaultTestLoader.loadTestsFromModule(testpolyhedra))
#         testsuite.addTests(unittest.defaultTestLoader.loadTestsFromModule(testph_types))
#         runner = unittest.TextTestRunner(verbosity=2)
#         runner.run(testsuite)
#     # testsuite.run(self)

# Execute all of the tests in all of the modules in the TgyCalc project
if __name__ == "__main__":
    localtime = time.asctime(time.localtime(time.time()))
    print("***** testall.py *****\n" + localtime + "************************")
    unittest.main(verbosity=2)

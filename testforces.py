import unittest
import forces
import math


class TestTowerTuneParams(unittest.TestCase):
    def testfunction(self):
        strut_count = 3
        level_count = 2
        h = [1.1, 1.2]
        # r = level_count * [[1, 1]]
        r = [[1.3, 1.4], [1.5, 1.6]]
        # l_twist = level_count * [math.pi * (1 / 2 - 1 / strut_count)]
        l_twist = [3.1, 3.2]
        # iface_twist = (level_count - 1) * [math.pi / strut_count]
        iface_twist = [4.1]
        iface_overlap = [0.3]
        t_params = forces.TowerTuneParams(n=strut_count, levels=level_count, height=h, radius=r, level_twist=l_twist,
                                          interface_twist=iface_twist, interface_overlap=iface_overlap)
        p_array = t_params.tune_param_array
        p_array += 0.1
        t_params.set_tune_params(p_array)
        new_p_array = t_params.tune_param_array
        for i, p in enumerate(p_array):
            self.assertEqual(p, new_p_array[i])
        strut_count = 3
        level_count = 3
        h = [1.1, 1.2, 1.3]
        # r = level_count * [[1, 1]]
        r = [[2.1, 2.2], [2.3, 2.4], [2.5, 2.6]]
        # l_twist = level_count * [math.pi * (1 / 2 - 1 / strut_count)]
        l_twist = [3.1, 3.2, 3.3]
        # iface_twist = (level_count - 1) * [math.pi / strut_count]
        iface_twist = [4.1, 4.2]
        iface_overlap = [0.3, 0.4]
        t_params = forces.TowerTuneParams(n=strut_count, levels=level_count, height=h, radius=r, level_twist=l_twist,
                                          interface_twist=iface_twist, interface_overlap=iface_overlap)
        p_array = t_params.tune_param_array
        p_array += 0.1
        t_params.set_tune_params(p_array)
        new_p_array = t_params.tune_param_array
        for i, p in enumerate(p_array):
            self.assertEqual(p, new_p_array[i])

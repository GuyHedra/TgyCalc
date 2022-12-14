import main
import unittest
# import testui


def prismtowerknownvalues():
    tower = main.prismTower(struts_per_level=3, level_count=1, overlap_ratio=0.1, heights=[1000], radii=[500, 500])
    print('strut length single level 3 strut', tower.struts[0].magnitude)
    assertEqual(1390, tower.struts[0].magnitude)
    assertEqual(1390, 5)


prismtowerknownvalues()

from unittest import TestCase

class TestprismTower(TestCase):
    def test_known_values(self):
        import main
        # 3 strut, one level
        prism = main.prismTower(struts_per_level=3,
                                level_count=1,
                                overlap_ratio=0.1,
                                heights=[1000],
                                radii=[500, 500],
                                hub_radius=12)
        self.assertEqual(1390, round(prism.struts[0].magnitude))
        self.assertEqual(1033, round(prism.tendon_verticals[0][0].magnitude))
        self.assertEqual(866, round(prism.tendon_waists[0][0].magnitude))
        self.assertEqual(866, round(prism.tendon_waists[1][0].magnitude))
        # 5 strut, one level
        prism = main.prismTower(struts_per_level=5,
                                level_count=1,
                                overlap_ratio=0.1,
                                heights=[1000],
                                radii=[500, 500],
                                hub_radius=12)
        self.assertEqual(1339, round(prism.struts[0].magnitude))
        self.assertEqual(1098, round(prism.tendon_verticals[0][0].magnitude))
        self.assertEqual(588, round(prism.tendon_waists[0][0].magnitude))
        self.assertEqual(588, round(prism.tendon_waists[1][0].magnitude))

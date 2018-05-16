# -*- coding: utf-8 -*-
"""
Created on 13 May 2018
@author: Dylan Jones

"""
from unittest import TestCase, main
from waveprop import Cell
from waveprop import BandFinder


class BandFinderTests(TestCase):

    def setUp(self):
        self.band_values = [14.79797, 26.46914]
        cell = Cell(10, 1, 0.8)
        self.bf = BandFinder(cell)

    def tearDown(self):
        self.bf = None

    def test_lower_limit(self):
        e = self.bf.lower_limit(0)
        self.assertAlmostEqual(e, self.band_values[0], 4)

    def test_upper_limit(self):
        e = self.bf.upper_limit(0)
        self.assertAlmostEqual(e, self.band_values[1], 4)


if __name__ == "__main__":
    main()

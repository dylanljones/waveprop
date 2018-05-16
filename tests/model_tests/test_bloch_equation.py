# -*- coding: utf-8 -*-
"""
Created on 13 May 2018
@author: Dylan

"""
from unittest import TestCase, main
from waveprop import Cell
from waveprop.model.bands import get_gaps, right_equation_from_e, right_equation, bloch_vector
import math


class BlochTests(TestCase):

    def setUp(self):
        v, a, d = 10, 1, 0.8
        self.e = 8
        self.cell = Cell(v, a, d)
        self.q = math.sqrt(2*self.e)
        self.kappa = math.sqrt(2*(v-self.e))

    def test_get_gaps(self):
        bands = [[1, 2], [4, 5]]
        gaps = get_gaps(bands)
        self.assertEqual(gaps, [[0, 1], [2, 4]])

    def test_right_equation(self):
        eq_r = right_equation(self.q, self.kappa, self.cell)
        self.assertAlmostEqual(eq_r, 0.51764+0j, 4)

    def test_right_equation_from_e(self):
        eq_r = right_equation_from_e(self.e, self.cell)
        self.assertAlmostEqual(eq_r, 0.51764+0j, 4)

    def test_bloch_vektor(self):
        k = bloch_vector(self.e, self.cell)
        self.assertAlmostEqual(k, 1.026705, 4)


if __name__ == "__main__":
    main()

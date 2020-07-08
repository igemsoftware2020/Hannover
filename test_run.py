#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a awesome
        python script!
"""

import numpy as np

# own modules
import biofilm


def test_evolution():
    initial_vals = {"position": np.array([4, 4]),  "length": 4,
                    "velocity": np.array([1, -1]), "angle": np.pi / 3, "index": 0}

    bfilm = biofilm.Biofilm()
    bfilm.initial_configuration(initial_vals)
    bfilm.print_current_state()
    for _ in range(0, 11):
        bfilm.evolve_bacteria()

    bfilm.plot_iterations()
    bfilm.print_current_state()


test_evolution()

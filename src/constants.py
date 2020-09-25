#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ********************************************************************************************
# imports
import os
from pathlib import Path

import numpy as np


class Constants:
    """
    This class is for managing and storing different biological and physical constants,
    which are used in the simulation
    """
    # FILE PATHS
    ROOT_DIRECTORY = Path(os.getcwd())
    OUTPUT_PATH = ROOT_DIRECTORY / 'output'
    SOURCE_PATH = ROOT_DIRECTORY / 'src'

    # SIMULATION PARAMETERS
    TIME_STEP = 1  # in seconds
    NUM_INITIAL_BACTERIA = 3

    # constants regarding B. subtilius
    BSUB_WIDTH = 1  # MICROMETERS https://en.wikipedia.org/wiki/Bacillus_subtilis
    BSUB_LENGTH = 4.9  # MICROMETERS https://en.wikipedia.org/wiki/Bacillus_subtilis
    BSUB_MASS = 10 ** (-15)  # kg
    BSUB_MORTALITY_RATE = 0.00
    BSUB_CRITICAL_LENGTH = 9  # MICROMETERS

    # SPEEDS
    BSUB_FREE_MEAN_SPEED = 50  # 50 MICROMETERS / s mean speed

    MAX_RADIAL_SPEED = 6  # [um / h] DOI 10.1126/science.abb8501 (2020).
    MAX_LATERAL_SPEED = 8  # [um / h] DOI 10.1126/science.abb8501 (2020)

    # ADHESION FORCES and VISCOSITY
    MAX_CELL_SUBSTRATE_ADHESION = 5.08 * 1E-9  # [N] DOI 10.1016/S0167-7012(99)00137-2
    MAX_CELL_CELL_ADHESION = 6.81 * 1E-9  # [N] DOI 10.1016/S0167-7012(99)00137-2
    EFFECTIVE_VISCOSITY_EPS = np.log(1E3)  # [Pa * s] : of bacterial P. aeruginosa PAO1 10.1103/PhysRevLett.93.098102

    # DOUBLING TIMES AND GROWTH RATES
    BSUB_DOUBLING_TIME = 720  # SECONDS DOI: 10.1128/jb.167.1.219-230.1986
    ECOLI_DOUBLING_TIME = 1200  # [s] LINK https://de.wikipedia.org/wiki/Generationszeit
    ECOLI_GROWTH_RATE = np.log(2) / ECOLI_DOUBLING_TIME  # [1 / s]
    BSUB_GROWTH_RATE = np.log(2) / BSUB_DOUBLING_TIME  # [1 / s]
    BSUB_GROWTH_FACTOR = 0.02  # [um / s]

    MOTION_ACTIVATION_PROBABILITY = 0.005
    MOTION_DEACTIVATION_PROBABILITY = 0.01

    # constants regarding visualisation
    WINDOW_SIZE = (500, 750)
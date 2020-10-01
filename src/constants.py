#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
# ********************************************************************************************
# imports
import os
from pathlib import Path
from tkinter import filedialog


class Constants:
    """
    This class is for managing and storing different biological and physical constants,
    which are used in the simulation
    """

    def __init__(self):
        # FILE PATHS
        self.root_dir = Path(os.getcwd())
        self.output_path = self.root_dir / 'output'
        self.source_path = self.root_dir / 'src'

        # SIMULATION PARAMETERS
        self.simulation_constants = {}
        self.duration = 10
        self.num_initial_bac = 0
        self.time_step = 1
        self.window_size = (500, 750)

        # Bacteria constants
        self.bac_type = None
        self.bac_constants = {}

    def set_bacteria_constants(self, default=True):
        if self.bac_type == "B.Sub." and default:
            self.bac_constants = Constants.get_bsub_constants()
        elif self.bac_type == "E. Coli." and default:
            self.bac_constants = Constants.get_ecoli_constants()

    def get_bacteria_constants(self):
        if self.bac_constant:
            return self.bac_constants
        else:
            raise ValueError("Bacteria constants are not set.")

    def get_simulation_constants(self):
        if self.simulation_constants:
            return self.simulation_constants
        else:
            raise ValueError(" constants are not set.")

    def get_root_dir(self):
        return self.root_dir

    def set_root_dir(self):
        new_root = filedialog.askdirectory(initialdir=self.get_root_dir())
        self.root_dir = os.chdir(new_root)
        self.output_path = self.root_dir / 'output'
        self.source_path = self.root_dir / 'src'

    def set_simulation_constants(self):
        self.simulation_constants = {
            "Num_initial_bacteria": self.num_initial_bac,
            "Window_size": self.window_size,
            "Time_step": self.time_step,
            "Duration": self.duration
        }

    @staticmethod
    def get_bsub_constants():
        return {
            "LENGTH": 4.9,  # [um] https://en.wikipedia.org/wiki/Bacillus_subtilis
            "WIDTH": 1,  # [um] https://en.wikipedia.org/wiki/Bacillus_subtilis
            "MASS": 10 ** (-15),  # [kg]
            "MORTALITY_RATE": 0.00,
            "CRITICAL_LENGTH": 9,  # [um]
            "FREE_MEAN_SPEED": 50,  # [um / s]
            "DOUBLING_TIME": 720,  # [s] DOI: 10.1128/jb.167.1.219-230.1986
            "GROWTH_RATE": np.log(2) / 720,  # [1 / s]
        }

    @staticmethod
    def get_ecoli_constants():
        return {
            "LENGTH": 2,  # [um] https://en.wikipedia.org/wiki/Escherichia_coli
            "WIDTH": 0.5,  # [um] https://en.wikipedia.org/wiki/Escherichia_coli
            "MASS": 10 ** (-15),  # [kg]
            "MORTALITY_RATE": 0.00,
            "CRITICAL_LENGTH": 9,  # [um]
            "FREE_MEAN_SPEED": 50,  # [um / s]
            "DOUBLING_TIME": 1200,  # [s] DOI: 10.1128/jb.167.1.219-230.1986
            "GROWTH_RATE": np.log(2) / 1200,  # [1 / s]
        }

    # SPEEDS
    BSUB_FREE_MEAN_SPEED = 50  # 50 MICROMETERS / s mean speed

    MAX_RADIAL_SPEED = 6  # [um / h] DOI 10.1126/science.abb8501 (2020).
    MAX_LATERAL_SPEED = 8  # [um / h] DOI 10.1126/science.abb8501 (2020)

    # ADHESION FORCES and VISCOSITY
    MAX_CELL_SUBSTRATE_ADHESION = 5.08 * 1E-9  # [N] DOI 10.1016/S0167-7012(99)00137-2
    MAX_CELL_CELL_ADHESION = 6.81 * 1E-9  # [N] DOI 10.1016/S0167-7012(99)00137-2
    EFFECTIVE_VISCOSITY_EPS = np.log(1E3)  # [Pa * s] : of bacterial P. aeruginosa PAO1 10.1103/PhysRevLett.93.098102

    MOTION_ACTIVATION_PROBABILITY = 0.005
    MOTION_DEACTIVATION_PROBABILITY = 0.01

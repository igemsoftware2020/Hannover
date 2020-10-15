#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
# ********************************************************************************************
# imports
import os
from datetime import datetime
from pathlib import Path
from tkinter import filedialog
from typing import Dict


class Constants:
    """
    This class is for managing and storing the different biological and physical constant,
    which are used in the simulation
    """

    # Global parameters
    # ADHESION FORCES and VISCOSITY
    MAX_CELL_SUBSTRATE_ADHESION = 7 * 1E-9  # 5.08 * 1E-9  # [N] DOI 10.1016/S0167-7012(99)00137-2
    MAX_CELL_CELL_ADHESION = 6.81 * 1E-9  # [N] DOI 10.1016/S0167-7012(99)00137-2
    EFFECTIVE_VISCOSITY_EPS = np.log(1E3)  # # [Pa * s] : of bacterial P. aeruginosa PAO1 10.1103/PhysRevLett.93.098102
    EFFECTIVE_VISCOSITY_H2O = 0.7805 * 1E-3  # [Pa * s]: at ~ 30 °C https://wiki.anton-paar.com/en/water/

    # MAX_RADIAL_SPEED = 6  # [um / h] DOI 10.1126/science.abb8501 (2020).
    # MAX_LATERAL_SPEED = 8  # [um / h] DOI 10.1126/science.abb8501 (2020)

    def __init__(self, bac_type: str):
        # FILE PATHS
        self.root_dir = ""
        self.output_path = ""
        self.info_path = ""

        # SIMULATION PARAMETERS
        self.num_initial_bac = 3
        self.time_step = 1
        self.window_size = (450, 900)
        self.duration = 60  # simulation time in minutes
        self.sim_dict = {}

        # Bacteria constant
        self.bac_type = bac_type
        self.bac_constants = {}

    def __repr__(self):
        bac_dict = self.get_bac_constants()
        sim_dict = self.get_simulation_constants()
        paths = self.get_paths()

        def append_dic_str(s: str, d: Dict):
            for key, values in d.items():
                s += f"{str(key)} :   {str(values)}\n"
            return s

        repr_str = f"\n ******  PATHS   ******\n "
        repr_str = append_dic_str(repr_str, paths)
        repr_str += f"\n ******  CONSTANTS   ******\n   (check documentation for units)\n\n "
        repr_str += f"* Constants of {self.bac_type} *\n"
        repr_str = append_dic_str(repr_str, bac_dict)
        repr_str += f"\n * Simulation constant *\n"
        repr_str = append_dic_str(repr_str, sim_dict)
        return repr_str

    def set_bacteria_constants(self, default=True):
        """ set constants according to selected bacteria strain """
        if self.bac_type == "B.Sub." and default:
            self.bac_constants = Constants.get_bsub_constants()
        elif self.bac_type == "E.Coli." and default:
            self.bac_constants = Constants.get_ecoli_constants()

    def get_bac_constants(self, key: str = None):
        """ return used constants as dict. If key is given, returns the respective constants."""
        dict_constants = self.bac_constants
        if key and (key in dict_constants.keys()):
            return dict_constants[key]
        else:
            return dict_constants

    def get_paths(self, key=None):
        """
        returns root, output and path of info file in a dictionary.
        If key is given, return respective path from dictionary
        """
        paths_dir = {
            "root": self.root_dir,
            "output": self.output_path,
            "info": self.info_path
        }
        if key:
            if key not in paths_dir.keys():
                return paths_dir
            elif key in paths_dir.keys():
                return paths_dir[key]
        else:
            return paths_dir

    def set_paths(self, default: bool = True):
        """
        Sets paths used for storing log file and plots.
        Default root path is the working directory
         """
        if not default:
            path = Path(filedialog.askdirectory())
            os.chdir(path)
        else:
            path = Path(os.getcwd())

        self.root_dir = path
        self.output_path = self.root_dir / 'output'

        date_time = str(datetime.now().day) + str(datetime.now().month) + str(datetime.now().year) \
                    + '_' + str(datetime.now().hour) + 'h' + str(datetime.now().minute) + 'min'

        path_out = self.output_path / f'log_{date_time}'
        if not self.output_path.exists():
            self.output_path.mkdir()
        elif not path_out.exists():
            path_out.mkdir()
        self.info_path = path_out / f'log_{date_time}.json'

    def get_simulation_constants(self, key: str = None):
        """
        returns simulation constants in a dictionary.
        If key is given, return respective path from dictionary
        """
        dict_constants = self.sim_dict
        if key and (key in dict_constants.keys()):
            return dict_constants[key]
        else:
            return dict_constants

    def set_simulation_constants(self):
        """
        sets the simulation constants fix.
        """
        sim_dict = {
            "num_initial": self.num_initial_bac,
            "time_step": self.time_step,
            "window_size": self.window_size,
            "duration": self.duration
        }
        self.sim_dict = sim_dict

    @staticmethod
    def get_bsub_constants(key: str = None):
        """
        returns constants regarding Bacillus subtilis strain
        If key is given, return respective path from dictionary
        """
        bsub_dic = {
            "LENGTH": np.random.normal(loc=2.5, scale=2.5 * 0.14),
            "WIDTH": 1,  # [um] https://en.wikipedia.org/wiki/Bacillus_subtilis
            "MASS": 10 ** (-12),  # [kg]
            "MORTALITY_RATE": 0.0,
            "CRITICAL_LENGTH": 4.7,  # [um]
            "FREE_MEAN_SPEED": 8 / (60 * 60),  # [um / s]
            "DOUBLING_TIME": 7200,  # [s] DOI: 10.1128/jb.167.1.219-230.1986
            "GROWTH_RATE": 2.2 / 7200,  # [um / s]
            "MOTION_ACTIVATION_PROBABILITY": 0.005,
            "MOTION_DEACTIVATION_PROBABILITY": 0.01
        }
        if key and (key in bsub_dic):
            return bsub_dic[key]
        else:
            return bsub_dic

    @staticmethod
    def get_ecoli_constants(key: str = None):
        """
         returns constants regarding Escherichia coli strain
         If key is given, return respective path from dictionary
        """
        ecoli_dic = {
            "LENGTH": np.random.normal(loc=1, scale=1 * 0.14),  # [um] https://en.wikipedia.org/wiki/Escherichia_coli
            "WIDTH": 0.5,  # [um] https://en.wikipedia.org/wiki/Escherichia_coli
            "MASS": 10 ** (-12),  # [kg]
            "MORTALITY_RATE": 0.01,
            "CRITICAL_LENGTH": 2,  # [um]
            "FREE_MEAN_SPEED": 50,  # [um / s]
            "DOUBLING_TIME": 1200,  # [s] DOI: 10.1128/jb.167.1.219-230.1986
            "GROWTH_RATE": 1 / 1200,  # [um / s]
            "MOTION_ACTIVATION_PROBABILITY": 0.005,
            "MOTION_DEACTIVATION_PROBABILITY": 0.01
        }
        if key and (key in ecoli_dic):
            return ecoli_dic[key]
        else:
            return ecoli_dic

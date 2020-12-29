#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ********************************************************************************************
# imports
import time
from typing import Dict
import numpy as np
# custom libraries
from BiofilmSimulation.constants import Constants


# ********************************************************************************************
# Utility functions

def simulation_duration(func):
    def inner1(*args, **kwargs):
        # storing time before function execution
        begin = time.time()
        func(*args, **kwargs)
        # storing time after function execution
        end = time.time()
        print(f'Duration of {func.__name__} : {end - begin} s')

    return inner1


def prompt_log_at_start(constants: Constants):
    """ Log printed in terminal at start """
    print(f" ************ BIOFILM MODELING ************ \n"
          " A project of the iGEM Teams Hannover x Darmstadt\n")
    print(constants)


def print_dic(dic: Dict):
    for key, item in dic.items():
        print(f"  {key} :  {item}")


def get_grid_coordinates(constants: Constants, distance: float = 1.):
    """
    Creates a grid with limits according to the window size defined in constants.
    Returns equally spaced, discrete cartesian coordinates.
    The spacing between two points can be defined in distance (expects microns as a unit).
I think of the vector entries as the midpoint coordinates of cubes, which discretize the space.
    :param constants: Instance of constants class.
    :param distance: distance between two points in the grid
    """
    x_limit, y_limit, z_limit = constants.window_size
    x = np.arange(distance, x_limit - distance, distance)
    y = np.arange(distance, y_limit - distance, distance)
    z = np.arange(distance, z_limit - distance, distance)
    return x, y, z

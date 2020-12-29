#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ********************************************************************************************
# imports
import numpy as np
# custom libraries
from BiofilmSimulation.constants import Constants

# ********************************************************************************************


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

    x = list(np.arange(distance, x_limit - distance, distance))
    y = list(np.arange(distance, y_limit - distance, distance))
    z = list(np.arange(distance, z_limit - distance, distance))

    x = [round(element, 1) for element in x]
    y = [round(element, 1) for element in y]
    z = [round(element, 1) for element in z]
    return x, y, z


def get_empty_grid(grid):
    x_empty = [None for _ in range(len(grid[0]))]
    y_empty = [None for _ in range(len(grid[1]))]
    z_empty = [None for _ in range(len(grid[2]))]
    return x_empty, y_empty, z_empty


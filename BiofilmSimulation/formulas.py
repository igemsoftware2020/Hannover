#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.spatial.transform import Rotation as R


# ********************************************************************************************
# In this script, the formulas for the forces and rotations are stored

def stokes_drag_force(radius: float, velocity: np.ndarray, viscosity: float) -> np.ndarray:
    # Calculates Stokes' drag for a sphere with Reynolds number < 1.
    # [um * Pa * s  * um / s] = [um * kg / (m * s ** 2) * s  * um / s]
    # changed units to N
    return - 6 * np.pi * radius * viscosity * 1E-12 * velocity


def gravitational_force(mass: float, distance_over_surface: float) -> np.ndarray:
    # calculates gravitational force on a mass
    # F = m * g * e_z
    # [kg * um / s ** 2]
    # changed units to N
    return mass * 9.81 * distance_over_surface * np.asarray([0, 0, -1])


def lennard_jones_force(r, f_min, r_min):
    """
    Calculates the lennard-jones force at r. Formula is derived by calculating the gradient of the
     (12, 6) lennard jones potential.
    :param r: value at which the function is evaluated
    :param f_min: value of the global minimum/maximum
    :param r_min: value at which the returned force is 0
    :return: force derived from the lennard- jones potential at value r
    """
    epsilon = f_min * (-169 * (r_min / (2 ** (1 / 6))) / (252 * (7 / 13) ** (1 / 6) * 2 ** (5 / 6)))
    sigma = r_min / (2 ** (1 / 6))
    return 48 * epsilon * np.power(sigma, 12) / np.power(r, 13) - 24 * epsilon * np.power(sigma, 6) / np.power(r, 7)


def get_euclid_norm(array):
    """ returns the norm of each vector in parameter array"""
    for i in range(len(array)):
        array[i] = np.linalg.norm(array[i])
    return array


def rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.
    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + np.cos(angle) * (px - ox) - np.sin(angle) * (py - oy)
    qy = oy + np.sin(angle) * (px - ox) + np.cos(angle) * (py - oy)
    return qx, qy


def rotation_matrix_x(theta: float):
    # return numpy array with rotation matrix around x axis with angle theta
    r = R.from_euler('x', theta, degrees=True)
    return r


def rotation_matrix_y(theta: float):
    # return numpy array with rotation matrix around y axis with angle theta
    r = R.from_euler('y', theta, degrees=True)
    return r


def rotation_matrix_z(theta: float):
    # return numpy array with rotation matrix around z axis with angle theta
    r = R.from_euler('z', theta, degrees=True)
    return r


def apply_rotation(vector: np.ndarray, matrix: R):
    return matrix.apply(vector)

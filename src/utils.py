#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ********************************************************************************************
# imports
import json
import time
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
# custom libraries
import src.constants as Constants
from scipy.spatial.transform import Rotation as R


def get_z(data):
    x = []
    y = []
    z = []
    for bac in data['position'].index:
        x.append(data['position'][bac][-1][0])
        y.append(data['position'][bac][-1][1])
        z.append(data['position'][bac][-1][2])
    a = np.vstack((x, y, z)).T
    pos = pd.DataFrame(a, columns=['x', 'y', 'z'])
    pos = pos.sort_values('z', ascending=False)
    pos = pos[0:int(np.round((len(pos) * 0.1)))].values
    x = pos[:, 0]
    y = pos[:, 1]
    z = pos[:, 2]
    return x, y, z


# ********************************************************************************************
# Data handling


def write_log_template(info_file_path, constants: Constants):
    """ saves a json template for saving bacteria parameters"""
    with open(info_file_path, 'w+') as json_file:
        data = {'BACTERIA': {}, 'CONSTANTS': {}}
        data['CONSTANTS'].update(constants.get_simulation_constants())
        json.dump(data, json_file, indent=4)


def read_in_log(info_file_path) -> Dict:
    """ returns data from a json file as dictionary"""
    with open(info_file_path, "r") as json_file:
        data = json.load(json_file)
    return data


def bacteria_as_pandas(info_file_path) -> pd.DataFrame:
    """ sorts bacteria parameters into DataFrame"""
    data = read_in_log(info_file_path)
    return pd.DataFrame(data['BACTERIA']).transpose()


def constants_as_pandas(info_file_path):
    """ sorts all saved constant into a Dataframe and returns"""
    data = read_in_log(info_file_path)
    return pd.DataFrame(data['CONSTANTS']).transpose()


def save_dict_as_json(data: Dict, info_file_path: Path):
    with open(str(info_file_path), 'w') as json_file:
        json.dump(data, json_file)


def get_data_to_parameter(data: pd.DataFrame, key: str, exact: bool = False):
    """
    Gets complete bacteria data as a DataFrame.
    Resorts data to parameters 'key' into another DataFrame:
        bacteria_0_parameter    ... bacteria_N_parameter
    0       100                         NaN
    1       12                          NaN
    2       13                          10
    .
    .
    .

    If the key is position or velocity, calculates the absolute value first.
    Returned Data is sorted according to the iteration step.
    """
    dic = {}
    for index, bacteria in data.iterrows():
        df_bac = pd.DataFrame(bacteria).loc[key, :]
        for vectors in df_bac:
            if (key == 'velocity' or key == 'position' or key == 'acceleration') & (exact is False):
                # calculate norm for velocity and position
                vectors = get_euclid_norm(vectors)
            dic.update({str(index) + '_' + key: vectors})
    df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in dic.items()]))

    def isnan(vector):
        if np.isnan(np.min(vector)) and exact:
            return False
        elif exact:
            return True
        elif not exact:
            return pd.isna(vector)

    df = df.transform(lambda x: sorted(x, key=isnan, reverse=True))
    return df


# ********************************************************************************************
# Formulas


def stokes_drag_force(radius: float, velocity: np.ndarray, viscosity: float) -> np.ndarray:
    # Calculates Stokes' drag for a sphere with Reynolds number < 1.
    # [um * Pa * s  * um / s] = [um * kg / (m * s ** 2) * s  * um / s]
    # changed units to N
    return - 6 * np.pi * radius * viscosity * 1E-12 * velocity


def gravitational_force(mass: float) -> np.ndarray:
    # calculates gravitational force on a mass
    # F = m * g * e_z
    # [kg * um / s ** 2]
    # changed units to N
    return mass * 9.81 * np.asarray([0, 0, -1])


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

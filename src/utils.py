import json
import time
from datetime import datetime
from pathlib import Path
from typing import Dict

import matplotlib.pyplot as plt
# ********************************************************************************************
# imports
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation as R

# custom libraries
from src.constants import Constants as C


def write_log_template(info_file_path):
    constants = C()
    with open(info_file_path, 'w+') as json_file:
        data = {'BACTERIA': {}, 'CONSTANTS': {}}
        members = [attr for attr in dir(constants) if
                   not callable(getattr(constants, attr)) and not attr.startswith("__")]
        constants_dic = {}
        for constant in members:
            constants_dic.update({constant: str(getattr(constants, constant))})
        data['CONSTANTS'].update(constants_dic)
        json.dump(data, json_file, indent=4)


def read_in_log(info_file_path) -> Dict:
    with open(info_file_path, "r") as json_file:
        data = json.load(json_file)
    return data


def bacteria_as_pandas(info_file_path) -> pd.DataFrame:
    data = read_in_log(info_file_path)
    return pd.DataFrame(data['BACTERIA']).transpose()


def constants_as_pandas(info_file_path):
    data = read_in_log(info_file_path)
    return pd.DataFrame(data['CONSTANTS']).transpose()


def save_dict_as_json(data: Dict, info_file_path: Path):
    with open(str(info_file_path), 'w') as json_file:
        json.dump(data, json_file)


def get_data_to_parameter(data: pd.DataFrame, key: str):
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
            if key == 'velocity' or key == 'position':
                # calculate norm for velocity and position
                vectors = get_euclid_norm(vectors)
            dic.update({str(index) + '_' + key: vectors})
    df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in dic.items()]))
    df = df.transform(lambda x: sorted(x, key=pd.isnull, reverse=True))
    return df


def get_euclid_norm(array):
    """ returns the norm of each vector in parameter array"""
    for i in range(len(array)):
        array[i] = np.linalg.norm(array[i])
    return array


def plot_velocities(data: pd.DataFrame, save_path: Path, save_fig: bool = False):
    """
    Plots velocities of each bacteria and the mean velocity of all bacteria
    over the iteration step.
    """
    plot_data = get_data_to_parameter(data, 'velocity')
    means = plot_data.mean(axis=1, skipna=True)
    fig, (ax1, ax2) = plt.subplots(1, 2)
    for bacteria in plot_data:
        ax1.plot(plot_data.loc[:, bacteria])

    ax1.set_title('Velocities')
    ax1.set_xlabel('Step')
    ax1.set_ylabel('velocity')
    ax2.plot(means)
    ax2.set_title('Mean Velocity')
    ax2.set_xlabel('Step')
    ax2.set_ylabel('velocity')
    plt.show()

    if save_fig:
        fig.savefig(save_path / 'velocity_plot.jpeg')


def plot_positions(data: pd.DataFrame, save_path: Path, save_fig: bool = False):
    """
    Plots positions (as lengths of location vectors) of each bacteria and the mean position of all bacteria
    over the iteration step.
    """
    plot_data = get_data_to_parameter(data, 'position')
    means = plot_data.mean(axis=1, skipna=True)
    fig, (ax1, ax2) = plt.subplots(1, 2)
    for bacteria in plot_data:
        ax1.plot(plot_data.loc[:, bacteria])

    ax1.set_title('Position')
    ax1.set_xlabel('Step')
    ax1.set_ylabel('distance')
    ax2.plot(means)
    ax2.set_title('Mean position')
    ax2.set_xlabel('Step')
    ax2.set_ylabel('distance')
    plt.show()

    if save_fig:
        fig.savefig(save_path / 'positions_plot.jpeg')


def plot_force(data: pd.DataFrame, save_path: Path, save_fig: bool = False):
    """
    Plots force acting on each bacteria and the mean force acting on all bacteria
    over the iteration step.
    """
    plot_data = get_data_to_parameter(data, 'total_force')
    means = plot_data.mean(axis=1, skipna=True)
    fig, (ax1, ax2) = plt.subplots(1, 2)
    for bacteria in plot_data:
        ax1.plot(plot_data.loc[:, bacteria])

    ax1.set_title('Total force')
    ax1.set_xlabel('Step')
    ax1.set_ylabel('force')
    ax2.plot(means)
    ax2.set_title('Mean force')
    ax2.set_xlabel('Step')
    ax2.set_ylabel('force')
    plt.show()

    if save_fig:
        fig.savefig(save_path / 'force_plot.jpeg')


def plot_size(data: pd.DataFrame, save_path: Path, save_fig: bool = False):
    """
    Plots force acting on each bacteria and the mean force acting on all bacteria
    over the iteration step.
    """
    width_data = get_data_to_parameter(data, 'width')
    length_data = get_data_to_parameter(data, 'length')
    width_means = width_data.mean(axis=1, skipna=True)
    length_means = length_data.mean(axis=1, skipna=True)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    for bacteria in width_data:
        ax1.plot(width_data.loc[:, bacteria])
        ax2.plot(length_data.loc[:, bacteria.replace('width', 'length')])
    ax1.set_title('width')
    ax1.set_xlabel('Step')
    ax1.set_ylabel('width')
    ax2.set_title('length')
    ax2.set_xlabel('Step')
    ax2.set_ylabel('length')
    ax3.plot(width_means)
    ax4.plot(length_means)
    ax3.set_title('width mean')
    ax3.set_xlabel('Step')
    ax3.set_ylabel('width')
    ax4.set_title('length mean')
    ax4.set_xlabel('Step')
    ax4.set_ylabel('length')
    plt.show()

    if save_fig:
        fig.savefig(save_path / 'size_plot.jpeg')


def get_info_file_path():
    date_time = str(datetime.now().hour) + 'h' + str(datetime.now().minute) + 'min_' + \
                str(datetime.now().day) + str(datetime.now().month) + \
                str(datetime.now().year)

    path_out = C.OUTPUT_PATH / f'log_{date_time}'
    path_out.mkdir()
    info_file_name = path_out / f'log_{date_time}.json'
    return info_file_name


def prompt_log_at_start(save_dir: str):
    return (f"********************* BIOFILM MODELING *********************\n"
            "NUMBER OF INITIAL BACTERIA\t {number_bacteria}\n"
            "NUMBER OF ITERATIONS\t {iterations}\n\n"
            "==================================================\n"
            "INITIAL DIMENSIONS (LENGTH, WIDTH)\t {BSUB_LENGTH},\t{BSUB_WIDTH}\n"
            "MASS\t {BSUB_MASS}\n"
            "GROWTH FACTOR\t {BSUB_GROWTH_FACTOR}\n"
            "CRITICAL LENGTH\t {BSUB_CRITICAL_LENGTH}\n\n"
            "SAVING AS \t {saving_dir}"
            .format(number_bacteria=C.START_NUMBER_BACTERIA, iterations=C.NUMBER_ITERATIONS,
                    type="B. subtilius", BSUB_LENGTH=C.BSUB_LENGTH,
                    BSUB_WIDTH=C.BSUB_WIDTH, BSUB_MASS=C.BSUB_MASS,
                    BSUB_CRITICAL_LENGTH=C.BSUB_CRITICAL_LENGTH,
                    BSUB_GROWTH_FACTOR=C.BSUB_GROWTH_FACTOR, saving_dir=save_dir))


def stokes_drag_force(radius: float, velocity: np.ndarray, viscosity=C.EFFECTIVE_VISCOSITY_EPS) -> np.ndarray:
    # Calculates Stokes' drag for a sphere with Reynolds number < 1.
    return - 6 * np.pi * radius * viscosity * velocity


def gravitational_force(mass: float) -> np.ndarray:
    # calculates gravitational force on a mass
    # F = m * g * e_z
    return mass * 9.81 * np.asarray([0, 0, -1])


def simulation_duration(func):
    def inner1(*args, **kwargs):
        # storing time before function execution
        begin = time.time()

        func(*args, **kwargs)

        # storing time after function execution
        end = time.time()
        print("Duration : ", func.__name__, end - begin)

    return inner1


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
    r = R.from_euler('x', theta)
    return r.as_matrix()


def rotation_matrix_y(theta: float):
    # return numpy array with rotation matrix around y axis with angle theta
    r = R.from_euler('y', theta)
    return r.as_matrix()


def rotation_matrix_z(theta: float):
    # return numpy array with rotation matrix around z axis with angle theta
    r = R.from_euler('z', theta)
    return r.as_matrix()


def apply_rotation(vector: np.ndarray, matrix: R):
    return matrix.apply(vector)

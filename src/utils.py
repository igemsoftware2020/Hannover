#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ********************************************************************************************
# imports
import json
import time
from pathlib import Path
from typing import Dict

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Ellipse
from scipy.spatial.transform import Rotation as R
from sklearn.linear_model import LinearRegression

# custom libraries
import src.constants as Constants


# ********************************************************************************************
# Functions for plotting data

def plot_velocities(data: pd.DataFrame, save_path: Path, save_fig: bool = False):
    """
    Plots velocities of each bacteria and the mean velocity of all bacteria
    over the iteration step.
    """
    plot_data = get_data_to_parameter(data, 'velocity')
    means = plot_data.mean(axis=1, skipna=True)
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.tight_layout()
    for bacteria in plot_data:
        ax1.plot(plot_data.loc[:, bacteria], '--', alpha=0.3)

    ax1.set_title('Velocities')
    ax1.set_xlabel('Time in s')
    ax1.set_ylabel('Velocity in um / s')
    ax2.plot(means)
    ax2.set_title('Mean Velocity')
    ax2.set_xlabel('Time in s')
    ax2.set_ylabel('Velocity in um / s')

    plt.ioff()

    if save_fig:
        path = Path(save_path).parent / 'velocity_plot.png'
        fig.savefig(str(path))
        plt.close(fig)
    else:
        plt.show()


def animate_positions(data: pd.DataFrame, save_path: Path, save_fig: bool = False):
    """
    Plots or saves (as mp4) an 2d animation of the biofilm.
    Animation is a top view of the biofilm and shows the trajectories of all bacteria in the simulation time.
    """
    plot_data = get_data_to_parameter(data, 'position', exact=True)
    living_data = get_data_to_parameter(data, 'living')

    fig, ax = plt.subplots()
    ax.set_xlabel("x / um")
    ax.set_ylabel("y / um")

    lines = []
    data = []
    living = []

    for bacteria in plot_data:
        x_data = [vector[0] for vector in plot_data[bacteria] if not np.isnan(np.min(vector))]
        y_data = [vector[1] for vector in plot_data[bacteria] if not np.isnan(np.min(vector))]
        lines.append(ax.plot(x_data, y_data, ), )
        data.append([x_data, y_data])
        living.append([living_data[bacteria.replace('position', 'living')]])

    def update(num, line_plots, dataLines, living_data):
        for line, dataLine, alive in zip(line_plots, dataLines, living_data):
            # update data for line plot: dataLine[0] = x data, dataLine[1] y data
            line[0].set_data(dataLine[0][num - 2:num], dataLine[1][num - 2:num])
            if alive[0] is False:
                line[0].set_color('black')
                line[0].set_alpha(0.8)
            ax.set_title(f"Trajectory of bacteria\npassed time: {round(num / 60, 2)} min")

        return lines,

    anim = animation.FuncAnimation(fig, update, frames=len(plot_data['bacteria_0_position']),

                                   interval=1000, repeat=False, fargs=[lines, data, living])

    plt.ioff()

    if save_fig:
        writer = animation.FFMpegWriter(fps=30, metadata=dict(artist='Me'), bitrate=-1)
        path = Path(save_path).parent / '2d_animation.mp4'
        anim.save(str(path), writer=writer)
        plt.close(fig)
    else:
        plt.show()


def animate_3d(data: pd.DataFrame, save_path: Path, save_fig: bool = False):
    """
    Plots or saves (as mp4) an 3d animation of the biofilm.
    Shows the trajectories of all bacteria in the simulation time in 3 dimensions.
    """
    np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
    plot_data = get_data_to_parameter(data, 'position', exact=True)

    fig = plt.figure()
    ax = p3.Axes3D(fig)
    ax.set_xlabel("x / um")
    ax.set_ylabel("y / um")
    ax.set_zlabel("z / um")

    lines = []
    data = []
    for bacteria in plot_data:
        x_data = np.asarray([vector[0] for vector in plot_data[bacteria] if not np.isnan(np.min(vector))])
        y_data = np.asarray([vector[1] for vector in plot_data[bacteria] if not np.isnan(np.min(vector))])
        z_data = np.asarray([vector[2] for vector in plot_data[bacteria] if not np.isnan(np.min(vector))])
        lines.append(ax.plot(x_data, y_data, z_data, alpha=0.9), )
        data.append([x_data, y_data, z_data])
    lines = np.asarray(lines)
    data = np.asarray(data)

    def update(num, line_plots, dataLines):
        for line, dataLine in zip(line_plots, dataLines):
            # update data for line plot: dataLine[0] = x data, dataLine[1] y data
            line[0].set_data(dataLine[0][num - 2:num], dataLine[1][num - 2:num])
            line[0].set_3d_properties(dataLine[2][num - 2:num])
            ax.set_title(f"Trajectory of bacteria\npassed time: {round(num / 60, 2)} min")
        return lines,

    anim = animation.FuncAnimation(fig, update, frames=len(plot_data['bacteria_0_position']),
                                   interval=100, repeat=False, fargs=[lines, data])

    plt.ioff()
    if save_fig:
        writer = animation.FFMpegWriter(fps=50, metadata=dict(artist='Me'), bitrate=-1)
        path = Path(save_path).parent / '3d_animation.mp4'
        anim.save(str(path), writer=writer)
        plt.close(fig)
    else:
        plt.show()


def plot_as_ellipse(data: pd.DataFrame, save_path: Path, save_fig: bool = False):
    pos_data = get_data_to_parameter(data, 'position')
    angle_data = get_data_to_parameter(data, 'angle')
    length_data = get_data_to_parameter(data, 'length')
    fig, ax = plt.subplots(1, 1)
    fig.tight_layout()
    ells = []
    for bacteria in pos_data:
        patch = Ellipse(xy=pos_data.loc[0, bacteria], width=1,
                        height=length_data.loc[0, bacteria],
                        angle=angle_data.loc[0, bacteria][0])
        ells.append(patch)

    for e in ells:
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)
        e.set_alpha(0.4)
    plt.show()

    if save_fig:
        fig.savefig(save_path / 'ellipse_plot.jpeg')


def plot_positions(data: pd.DataFrame, save_path: Path, save_fig: bool = False):
    """
    Plots positions (as lengths of location vectors) of each bacteria and the distance over the surface.
    """
    position_data = get_data_to_parameter(data, 'position')
    position_means = position_data.mean(axis=1, skipna=True)
    height_data = get_data_to_parameter(data, 'height')
    height_means = height_data.mean(axis=1, skipna=True)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, constrained_layout=True)

    for bacteria in position_data:
        ax1.plot(position_data.loc[:, bacteria], '--', alpha=0.3)
        ax2.plot(height_data.loc[:, bacteria.replace('position', 'height')], '--', alpha=0.3)

    ax1.set_title('Distance from origin')
    ax1.set_xlabel('Time in s')
    ax1.set_ylabel('Distance in um')

    ax2.set_title('Distance from surface')
    ax2.set_xlabel('Time in s')
    ax2.set_ylabel('Height in um')

    ax3.plot(position_means)
    ax3.set_xlabel('Time in s')
    ax3.set_ylabel('Mean distance in um')

    ax4.plot(height_means)
    ax4.set_xlabel('Time in s')
    ax4.set_ylabel('Mean height in um')

    plt.ioff()
    if save_fig:
        path = Path(save_path.parent) / 'positions_plot.png'
        fig.savefig(path)
        plt.close(fig)
    else:
        plt.show()


def plot_force(data: pd.DataFrame, save_path: Path, save_fig: bool = False):
    """
    Plots force acting on each bacteria and the mean force acting on all bacteria
    over the iteration step. Also plots accelerations.
    """
    plot_data = get_data_to_parameter(data, 'total_force')
    acc_data = get_data_to_parameter(data, 'acceleration')
    force_mean = plot_data.mean(axis=1, skipna=True)
    acc_mean = acc_data.mean(axis=1, skipna=True)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, constrained_layout=True)

    for bacteria in plot_data:
        ax1.plot(plot_data.loc[:, bacteria], '.', alpha=0.3)

    ax1.set_title('Total force')
    ax1.set_xlabel('Time in s')
    ax1.set_ylabel('Force in N')

    ax2.plot(acc_data)
    ax2.set_title('Total acceleration')
    ax2.set_xlabel('Time in s')
    ax2.set_ylabel('Acceleration in um/s²')

    ax3.plot(force_mean, '.')
    ax3.set_xlabel('Time in s')
    ax3.set_ylabel('mean force in N')

    ax4.plot(acc_mean, '-')
    ax4.set_xlabel('Time in s')
    ax4.set_ylabel('mean acceleration in um/s²')

    plt.ioff()
    if save_fig:
        path = Path(save_path).parent / 'force_plot.png'
        fig.savefig(path)
        plt.close(fig)
    else:
        plt.show()


def plot_sizes(data: pd.DataFrame, save_path: Path, save_fig: bool = False):
    """
    Plots force acting on each bacteria and the mean force acting on all bacteria
    over the iteration step.
    """
    mass_data = get_data_to_parameter(data, 'mass')
    mass_mean = mass_data.mean(axis=1, skipna=True)
    length_data = get_data_to_parameter(data, 'length')
    length_means = length_data.mean(axis=1, skipna=True)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, constrained_layout=True)

    for bacteria in mass_data:
        ax1.plot(mass_data.loc[:, bacteria], '.', alpha=0.3)
        ax2.plot(length_data.loc[:, bacteria.replace('mass', 'length')], '.', alpha=0.3)
    ax1.set_title('Masses')
    ax1.set_xlabel('Time in s')
    ax1.set_ylabel('Mass in kg')

    ax2.set_title('Lengths')
    ax2.set_xlabel('Time in s')
    ax2.set_ylabel('length in um')

    ax3.plot(mass_mean, '.')
    ax3.set_title('Mean of bacteria masses')
    ax3.set_xlabel('Time in s')
    ax3.set_ylabel('Mass in kg')

    ax4.plot(length_means, '.')
    ax4.set_title('Mean of bacteria lengths')
    ax4.set_xlabel('Time in s')
    ax4.set_ylabel('mean length in um')

    plt.ioff()
    if save_fig:
        path = Path(save_path).parent / 'size_plot.png'
        fig.savefig(path)
        plt.close(fig)
    else:
        plt.show()


def plot_num(data: pd.DataFrame, save_path: Path, save_fig: bool = False):
    live = get_data_to_parameter(data, 'living')
    num = live[live == True].count(axis=1)
    x, y_fit, slope, generation_time = get_gent(data)
    '''plot data'''
    fig, (ax1, ax2) = plt.subplots(2, 1)
    fig.tight_layout()
    ax1.plot(num, color='b')
    ax1.set(xlabel='Time in s', ylabel='Bacteria Number', title='Bacteria Growth')
    ax2.plot(num, label='log curve')
    ax2.set(xlabel='Time in s', ylabel='Bacteria Number [log]', title='Bacteria Growth')
    ax2.plot(x, np.exp(y_fit), label='fit curve')
    ax2.legend(loc='lower right')
    ax2.text(0.1, 0.9, 'slope: ' + str(round(slope, 5)), transform=ax2.transAxes)
    ax2.text(0.1, 0.8, 'generation time: ' + str(round(generation_time, 5)), transform=ax2.transAxes)
    ax2.set_yscale('log')

    plt.ioff()
    if save_fig:
        path = Path(save_path).parent / 'growth_plot.png'
        plt.savefig(path)
        plt.close(fig)
    else:
        plt.show()


def dens_map(data: pd.DataFrame, save_path: Path, save_fig: bool = False):
    """Scatters the last positions of the bacteria and plots the density of bacteria. """
    x, y, z = last_pos(data)
    fig, (ax1, ax2) = plt.subplots(1, 2, constrained_layout=True)
    ax1.scatter(x, y, c='g', s=20, alpha=0.8, marker='x')
    sns.kdeplot(data=x, data2=y, ax=ax2, shade=True, cbar=False, cmap='mako', levels=200, thresh=0)
    plt.ioff()
    if save_fig:
        path = Path(save_path).parent / 'density_plot.png'
        plt.savefig(path)
        plt.close(fig)
    else:
        plt.show()


def scatter_last_positions(data: pd.DataFrame, save_path: Path, save_fig: bool = False):
    x, y, z = last_pos(data)
    X, Y, Z = get_z(data)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, c='g', s=50, alpha=0.8, marker='o')
    ax.plot_trisurf(X, Y, Z, cmap='Greens', edgecolor='none')
    ax.set_xlabel('x / um')
    ax.set_ylabel('y / um')
    ax.set_zlabel('z / um')
    ax.view_init(60)

    plt.ioff()
    if save_fig:
        path = Path(save_path).parent / 'scatter_last_plot.png'
        plt.savefig(path)
        plt.close(fig)
    else:
        plt.show()


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
    # [um * Pa * s 1/1E-6 * um / s] = [um * kg / (um * s ** 2) * s  * um / s] = [um kg / (s ** 2)]
    return - 6 * np.pi * radius * viscosity * 1E-12 * velocity


def gravitational_force(mass: float) -> np.ndarray:
    # calculates gravitational force on a mass
    # F = m * g * e_z
    # [kg * um / s ** 2]
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


def get_gent(data: pd.DataFrame):
    live = get_data_to_parameter(data, 'living')  # get data
    y = live[live == True].count(axis=1).values  # transform data and return an array
    y = y[y != y[0]]  # cut out bacteria in lag phase
    y = [np.log(value) if value > 0 else value for value in y]  # transform data

    x = live.index[
        live[live == True].count(axis=1) != live[live == True].count(axis=1)[0]].to_numpy()  # get index array
    '''start linear regression'''
    model = LinearRegression(fit_intercept=True)
    model.fit(x[:, np.newaxis], y)  # fit the data
    generation_time = np.log(2) / model.coef_[0]  # compute generation time
    y_fit = model.predict(x[:, np.newaxis])  # get fitted curve
    slope = model.coef_[0]  # slope is growth coefficient
    return x, y_fit, slope, generation_time


def last_pos(data):
    last_cord_x = []
    last_cord_y = []
    last_cord_z = []
    for bac in data['position'].index:
        last_cord_x.append(data['position'][bac][-1][0])
        last_cord_y.append(data['position'][bac][-1][1])
        last_cord_z.append(data['position'][bac][-1][2])
    return last_cord_x, last_cord_y, last_cord_z


def prompt_log_at_start(constants: Constants):
    """ Log printed in terminal at start """
    print(f" ************ BIOFILM MODELING ************ \n"
          " A project of the iGEM Teams Hannover x Darmstadt\n")
    print(constants)


def print_dic(dic: Dict):
    for key, item in dic.items():
        print(f"  {key} :  {item}")

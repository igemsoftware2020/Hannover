#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pathlib import Path
from BiofilmSimulation.constants import Constants
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import pandas as pd
import seaborn as sns
# custom libraries
from BiofilmSimulation.data_handling import get_data_to_parameter, get_z, constants_as_pandas
from sklearn.linear_model import LinearRegression


# ********************************************************************************************
# Functions for plotting data



def animate_positions(data: pd.DataFrame, save_path: Path, save_fig: bool = False, time_step: int = 1):
    """
    Plots or saves (as mp4) an 2d animation of the biofilm.
    Animation is a top view of the biofilm and shows the trajectories of all bacteria in the simulation time.
    """
    plot_data = get_data_to_parameter(data, 'position', exact=True)
    living_data = get_data_to_parameter(data, 'living')

    fig, ax = plt.subplots()
    ax.set_xlabel("x / um")
    ax.set_ylabel("y / um")
    ax.set_xlim(1000)
    ax.set_ylim(1000)

    lines = []
    data = []
    living = []

    path = Path(save_path).parent / Path('2d_animation.mp4')

    for bacteria in plot_data:
        x_data = np.asarray([vector[0] for vector in plot_data[bacteria] if not np.isnan(np.min(vector))])
        y_data = np.asarray([vector[1] for vector in plot_data[bacteria] if not np.isnan(np.min(vector))])
        lines.append(ax.plot(x_data, y_data, ), )
        data.append([x_data, y_data])
        living.append([living_data[bacteria.replace('position', 'living')]])

    def update(num, line_plots, data_lines, data_living):
        for line, dataLine, alive in zip(line_plots, data_lines, data_living):
            # update data for line plot: dataLine[0] = x data, dataLine[1] y data
            line[0].set_data(dataLine[0][num - 5:num], dataLine[1][num - 5:num])
            if alive[0] is False:
                line[0].set_color('black')
                line[0].set_alpha(0.8)
            ax.set_title(f"Trajectory of bacteria\npassed time: - {round(num * time_step / 60, 2)} min")
        return lines,

    anim = animation.FuncAnimation(fig, update, frames=len(plot_data['bacteria_0_position']),
                                   interval=500, repeat=False, fargs=[lines, data, living])

    plt.ioff()

    if save_fig:
        writer = animation.FFMpegWriter(fps=30, metadata=dict(artist='Me'), bitrate=-1)
        anim.save(path, writer=writer)
        plt.close(fig)
    else:
        plt.show()


def animate_3d(data: pd.DataFrame, save_path: Path, save_fig: bool = False, time_step: int = 1):
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
    ax.set_xlim(1000)
    ax.set_ylim(1000)
    ax.set_zlim(50)
    lines = []
    data = []
    for bacteria in plot_data:
        x_data = (np.asarray([vector[0] for vector in plot_data[bacteria] if not np.isnan(np.min(vector))]))
        y_data = (np.asarray([vector[1] for vector in plot_data[bacteria] if not np.isnan(np.min(vector))]))
        z_data = (np.asarray([vector[2] for vector in plot_data[bacteria] if not np.isnan(np.min(vector))]))
        lines.append(ax.plot(x_data, y_data, z_data, alpha=0.8), )
        data.append([x_data, y_data, z_data])
    lines = np.asarray(lines)
    data = np.asarray(data)

    def update(num, line_plots, dataLines):
        for line, dataLine in zip(line_plots, dataLines):
            # update data for line plot: dataLine[0] = x data, dataLine[1] y data
            line[0].set_data(dataLine[0][num - 5:num], dataLine[1][num - 5:num])
            line[0].set_3d_properties(dataLine[2][num - 5:num])
            ax.set_title(f"Trajectory of bacteria\npassed time: {round(num * time_step/ 60, 2)} min\n"
                         )
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


def plot_velocities(data: pd.DataFrame, save_path: Path, save_fig: bool = False, time_step: int = 1):
    """
    Plots velocities of each bacteria and the mean velocity of all bacteria
    over the iteration step.
    """
    plot_data = get_data_to_parameter(data, 'velocity')
    means = plot_data.mean(axis=1, skipna=True)
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.tight_layout()
    for bacteria in plot_data:
        ax1.plot(plot_data.loc[:, bacteria].index * time_step,
                 plot_data.loc[:, bacteria].values,
                 '--', alpha=0.3
                 )

    ax1.set_title('Velocities')
    ax1.set_xlabel('Time in s')
    ax1.set_ylabel('Velocity in um / s')

    ax2.plot(means.index * time_step, means.values)
    ax2.set_title('Mean Velocity')
    ax2.set_xlabel('Time in s')
    ax2.set_ylabel('Velocity in um / s')

    plt.ioff()

    if save_fig and save_path:
        path = Path(save_path).parent / 'velocity_plot.png'
        fig.savefig(str(path))
        plt.close(fig)
    else:
        plt.show()


def plot_positions(data: pd.DataFrame, save_path: Path, save_fig: bool = False, time_step: int = 1):
    """
    Plots positions (as lengths of location vectors) of each bacteria and the distance over the surface.
    """
    position_data = get_data_to_parameter(data, 'position')
    position_means = position_data.mean(axis=1, skipna=True)
    height_data = get_data_to_parameter(data, 'height')
    height_means = height_data.mean(axis=1, skipna=True)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, constrained_layout=True)

    for bacteria in position_data:
        ax1.plot(position_data.loc[:, bacteria].index * time_step,
                 position_data.loc[:, bacteria].values,
                 '--', alpha=0.3
                 )
        ax2.plot(height_data.loc[:, bacteria.replace('position', 'height')].index * time_step,
                 height_data.loc[:, bacteria.replace('position', 'height')].values,
                 '--', alpha=0.3
                 )

    ax1.set_title('Distance from origin')
    ax1.set_xlabel('Time in s')
    ax1.set_ylabel('Distance in um')

    ax2.set_title('Distance from surface')
    ax2.set_xlabel('Time in s')
    ax2.set_ylabel('Height in um')

    ax3.plot(position_means.index * time_step, position_means)
    ax3.set_xlabel('Time in s')
    ax3.set_ylabel('Mean distance in um')

    ax4.plot(height_means.index * time_step, height_means)
    ax4.set_xlabel('Time in s')
    ax4.set_ylabel('Mean height in um')

    plt.ioff()
    if save_fig:
        path = Path(save_path).parent / 'positions_plot.png'
        fig.savefig(path)
        plt.close(fig)
    else:
        plt.show()


def plot_force(data: pd.DataFrame, save_path: Path, save_fig: bool = False, time_step: int = 1):
    """
    Plots force acting on each bacteria and the mean force acting on all bacteria
    over the iteration step. Also plots accelerations.
    """

    plot_data = get_data_to_parameter(data, 'total_force') * 1e9
    acc_data = get_data_to_parameter(data, 'acceleration')
    force_mean = plot_data.mean(axis=1, skipna=True)
    acc_mean = acc_data.mean(axis=1, skipna=True)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, constrained_layout=True)

    for bacteria in plot_data:
        ax1.plot(plot_data.loc[:, bacteria].index * time_step,
                 plot_data.loc[:, bacteria].values
                 , '.', alpha=0.3, markersize=3
                 )

    ax1.set_title('Total force')
    ax1.set_xlabel('Time in s')
    ax1.set_ylabel('Force in nN')

    ax2.plot(acc_data.index * time_step, acc_data.values, '.', alpha=0.3, markersize=1)
    ax2.set_title('Total acceleration')
    ax2.set_xlabel('Time in s')
    ax2.set_ylabel('Acceleration in um/s²')

    ax3.plot(force_mean.index * time_step, force_mean.values, '-')
    ax3.set_xlabel('Time in s')
    ax3.set_ylabel('mean force in nN')

    ax4.plot(acc_mean.index * time_step, acc_mean.values, '-')
    ax4.set_xlabel('Time in s')
    ax4.set_ylabel('mean acceleration in um/s²')

    plt.ioff()
    if save_fig:
        path = Path(save_path).parent / 'force_plot.png'
        fig.savefig(path)
        plt.close(fig)
    else:
        plt.show()


def plot_sizes(data: pd.DataFrame, save_path: Path, save_fig: bool = False, time_step: int = 1):
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
        ax1.plot(mass_data.loc[:, bacteria].index * time_step,
                 mass_data.loc[:, bacteria].values,
                 '.', alpha=0.3, markersize=1
                 )
        ax2.plot(length_data.loc[:, bacteria.replace('mass', 'length')].index * time_step,
                 length_data.loc[:, bacteria.replace('mass', 'length')].values,
                 '.', alpha=0.3, markersize=1
                 )

    ax1.set_title('Masses')
    ax1.set_xlabel('Time in s')
    ax1.set_ylabel('Mass in kg')

    ax2.set_title('Lengths')
    ax2.set_xlabel('Time in s')
    ax2.set_ylabel('length in um')

    ax3.plot(mass_mean.index * time_step, mass_mean.values, '.', markersize=3)
    ax3.set_title('Mean of bacteria masses')
    ax3.set_xlabel('Time in s')
    ax3.set_ylabel('Mass in kg')

    ax4.plot(length_means.index * time_step, length_means.values, '.', markersize=3)
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


def plot_num(data: pd.DataFrame, save_path: Path,  save_fig: bool = False, time_step: int = 1):

    live = get_data_to_parameter(data, 'living')
    num = live[live == True].count(axis=1)
    if get_gent(data, time_step) is None:
        return
    else:
        x, y_fit, slope, generation_time = get_gent(data, time_step)

    '''plot data'''
    fig, (ax1, ax2) = plt.subplots(2, 1)
    fig.tight_layout()
    ax1.plot(num.index * time_step, num.values, color='b')
    ax1.set(xlabel='Time in s', ylabel='Bacteria Number', title='Bacteria Growth')
    ax2.plot(num.index * time_step, num.values, label='log curve')
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


def get_gent(data: pd.DataFrame, time_step: int):

    live = get_data_to_parameter(data, 'living')  # get data

    y = live[live == True].count(axis=1).values  # transform data and return an array
    y = y[y != y[0]]  # cut out bacteria in lag phase
    y = [np.log(value) if value > 0 else value for value in y]  # transform data
    x = live.index[
        live[live == True].count(axis=1) != live[live == True].count(axis=1)[0]].to_numpy() * time_step # get index array
    '''start linear regression'''
    model = LinearRegression(fit_intercept=True)
    try:
        model.fit(x[:, np.newaxis], y)  # fit the data
    except ValueError:
        print("Simulation time to short to fit a regression on the bacteria number!")
        return None
    generation_time = np.log(2) / model.coef_[0]  # compute generation time
    y_fit = model.predict(x[:, np.newaxis])  # get fitted curve
    slope = model.coef_[0]  # slope is growth coefficient
    return x, y_fit, slope, generation_time


def last_pos(data):
    """ returns last coordinates of each bacteria in data"""
    last_cord_x = []
    last_cord_y = []
    last_cord_z = []
    for bac in data['position'].index:
        last_cord_x.append(data['position'][bac][-1][0])
        last_cord_y.append(data['position'][bac][-1][1])
        last_cord_z.append(data['position'][bac][-1][2])
    return last_cord_x, last_cord_y, last_cord_z


# Histogram PLots
def histo_length(data: pd.DataFrame, save_path: Path, save_fig: bool = False):
    data = data['length']
    length = []
    for bac in data.index:
        for i in data[bac]:
            length.append(i)
    end_length = np.array(length)  # the two for loops convert the dataFrame into a one dimensional array

    fig = sns.displot(end_length, kde=True,
                      **{'binwidth': 0.25, 'stat': 'density'})  # the histogram with a density function
    plt.xlabel('Bacteria length [um]')  # x label
    plt.ylabel('Normalized proportion')
    plt.tight_layout()
    plt.ioff()
    if save_fig:
        path = Path(save_path).parent / 'histo_length.png'
        plt.savefig(path)
    else:
        plt.show()


def histo_velocity(data: pd.DataFrame, save_path: Path, save_fig: bool = False):
    data = get_data_to_parameter(data, 'velocity')
    velocity = data.T.iloc[:, :].values.reshape(data.T.size)  # get the data as a one dimensional array

    end_velocity = velocity[np.logical_not(np.isnan(velocity))]  # filter all NaN values

    fig = sns.displot(end_velocity, kde=True,
                      **{'binwidth': 0.01, 'stat': 'density'})  # the histogram with a density function
    plt.xlabel('Bacteria velocity')
    plt.ylabel('Normalized proportion')

    plt.tight_layout()
    plt.ioff()
    if save_fig:
        path = Path(save_path).parent / 'histo_velocity.png'
        plt.savefig(path)
    else:
        plt.show()


def histo_force(data: pd.DataFrame, save_path: Path, save_fig: bool = False):
    plot_data = get_data_to_parameter(data, 'total_force') * 1e9

    force = plot_data.T.iloc[:, :].values.reshape(plot_data.T.size)  # get the data as a one dimensional array

    end_force = force[np.logical_not(np.isnan(force))]  # filter all NaN values

    fig = sns.displot(end_force, kde=True,
                      **{'binwidth': 0.1, 'stat': 'density'})  # the histogram with a density function
    plt.xlabel('Bacteria force')
    plt.tight_layout()
    plt.ylabel('Normalized proportion')

    plt.ioff()
    if save_fig:
        path = Path(save_path).parent / 'histo_force.png'
        plt.savefig(path)
    else:
        plt.show()


def lennard_jones_force_plot(r_min, f_min):
    """
    Forces resulting from the lennard jones potential. Reparameterized with r_min with F_LJ(r_min ) = 0
    and the absolute value of the global minimum f_min
    """

    def ljp(r, f_min, r_min):
        epsilon = f_min * (-169 * (r_min / (2 ** (1 / 6))) / (252 * (7 / 13) ** (1 / 6) * 2 ** (5 / 6)))
        sigma = r_min / (2 ** (1 / 6))
        return 48 * epsilon * np.power(sigma, 12) / np.power(r, 13) \
               - 24 * epsilon * np.power(sigma, 6) / np.power(r, 7)

    r = np.linspace(0.9999, 2, 1000)
    plt.plot(r, ljp(r, f_min, r_min))
    plt.xlabel('Distance in µm')
    plt.ylabel('Force in nN')
    plt.show()

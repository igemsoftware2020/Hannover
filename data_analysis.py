import os

from BiofilmSimulation.constants import Constants
from BiofilmSimulation.data_handling import bacteria_as_pandas, read_in_log, ask_for_log_dir
from BiofilmSimulation.plotting import histo_length, histo_velocity, histo_force
from BiofilmSimulation.plotting import plot_sizes, plot_force, plot_velocities, plot_positions, \
    plot_num, dens_map, animate_positions, animate_3d, lennard_jones_force_plot, plot_convex_hull, \
    scatter_last_positions, plot_delauny_triangulation, plot_surface
from pathlib import Path


def plotting(info_file_path):
    """ reads in data from info_file_path and plots data """
    data = bacteria_as_pandas(info_file_path)
    constants_log_path = (str(info_file_path).replace(".json", "_Constants.json"))
    constants = read_in_log(constants_log_path)
    time_step = constants['time_step']

    # Plot histograms
    histo_length(data, info_file_path, save_fig=True)
    histo_velocity(data, info_file_path, save_fig=True)
    histo_force(data, info_file_path, save_fig=True)

    # Distribution of biofilm on surface
    dens_map(data, info_file_path, save_fig=True)

    # Time series plots
    plot_num(data, info_file_path, time_step=time_step, save_fig=True)
    plot_velocities(data, info_file_path, time_step=time_step, save_fig=True)
    plot_positions(data, info_file_path, time_step=time_step, save_fig=True)
    plot_force(data, info_file_path, time_step=time_step, save_fig=True)
    plot_sizes(data, info_file_path, time_step=time_step, save_fig=True)

    # Animations
    data = bacteria_as_pandas(info_file_path)
    animate_positions(data, info_file_path, time_step=time_step, save_fig=True)


if __name__ == "__main__":
    save_fp = ask_for_log_dir()
    data = read_in_log(save_fp)
    plotting(save_fp)

#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ********************************************************************************************
# custom libraries
from src.biofilm import Biofilm
from src.constants import Constants
from src.data_handling import bacteria_as_pandas
from src.plotting import plot_sizes, plot_force, plot_velocities, plot_positions, animate_positions, plot_num, dens_map
from src.utils import prompt_log_at_start


def start_run(constant: Constants):
    """
    Starts a simulation run with the constants specified in constant
    :return info_file_path: path of log file
    """
    # Init a new Biofilm with above constant
    biofilm = Biofilm()
    # pass constant to biofilm object
    biofilm.constants = constant

    # Logging at begin
    prompt_log_at_start(biofilm.constants)
    # Save log file for
    info_file_path = constant.get_paths(key="info")

    info_file_path = constants.get_paths(key="info")
    biofilm.simulate_multiprocessing()
    plotting(info_file_path)


def plotting(info_file_path):
    """ reads in data from info_file_path and plots data """
    data = bacteria_as_pandas(info_file_path)
    plot_num(data, info_file_path, save_fig=True)
    dens_map(data, info_file_path, save_fig=True)
    plot_velocities(data, info_file_path, save_fig=True)
    plot_positions(data, info_file_path, save_fig=True)
    plot_force(data, info_file_path, save_fig=True)
    plot_sizes(data, info_file_path, save_fig=True)
    data = bacteria_as_pandas(info_file_path)
    animate_positions(data, info_file_path, save_fig=True)
    # animate_3d(data, info_file_path, save_fig=False)
    
# ********************************************************************************************
# main-method to start the program
# ********************************************************************************************


if __name__ == "__main__":
    # Set constant for modelling run
    constants = Constants(bac_type="B.Sub.")
    constants.num_initial_bac = 100
    constants.duration = 10
    constants.time_step = 0.5
    constants.set_paths(default=False)
    constants.window_size = (2000, 2000)
    constants.set_bacteria_constants()
    constants.set_simulation_constants()
    constants.set_paths(default=True)

    start_run(constants)

    # path = ask_for_log_dir()
    # plotting(path)

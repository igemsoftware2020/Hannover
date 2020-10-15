#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ********************************************************************************************
# imports

from pathlib import Path

# custom libraries
from src.biofilm import Biofilm
from src.constants import Constants
from src.utils import plot_sizes, plot_force, plot_velocities, plot_positions, bacteria_as_pandas, \
    prompt_log_at_start, animate_positions, animate_3d, plot_num, dens_map


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
    #biofilm.simulate()
    biofilm.simulate_multiprocessing()
    return info_file_path


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
    animate_3d(data, info_file_path, save_fig=False)
    
# ********************************************************************************************
# main-method to start the program
# ********************************************************************************************


if __name__ == "__main__":
    # Set constant for modelling run
    constants = Constants(bac_type="B.Sub.")
    constants.num_initial_bac = 69
    constants.duration = 6
    constants.time_step = 1
    constants.window_size = (2000, 2000)
    constants.set_bacteria_constants()
    constants.set_simulation_constants()
    constants.set_paths(default=True)

    #path = start_run(constants)
    path = '/home/david/IntelliJ_projects/iGEM-biofilm-model/output/log_16102020_1h34min/log_16102020_1h34min.json'
    plotting(Path(path))

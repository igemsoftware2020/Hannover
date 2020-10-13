#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ********************************************************************************************
# imports

from pathlib import Path

# custom libraries
from src.biofilm import Biofilm
from src.constants import Constants
from src.utils import plot_size, plot_force, plot_velocities, plot_positions, bacteria_as_pandas, \
    prompt_log_at_start, plot_num, dens_map, animate_positions


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
    biofilm.simulate(save_name=info_file_path)

    return info_file_path


def plotting(info_file_path):
    data = bacteria_as_pandas(info_file_path)
    plot_num(data, info_file_path, save_fig=True)
    dens_map(data, info_file_path, save_fig=True)
    plot_velocities(data, info_file_path, save_fig=False)
    plot_positions(data, info_file_path, save_fig=False)  # this one messes with data
    plot_force(data, info_file_path, save_fig=False)
    plot_size(data, info_file_path, save_fig=False)
    data = bacteria_as_pandas(info_file_path)
    animate_positions(data, info_file_path, save_fig=True)
    # animate_3d(data, info_file_path, save_fig=True)
    
# ********************************************************************************************
# main-method to start the program
# ********************************************************************************************


if __name__ == "__main__":
    # Set constant for modelling run
    constants = Constants(bac_type="B.Sub.")
    constants.num_initial_bac = 10
    constants.duration = 15
    constants.window_size = (2000, 2000)
    constants.set_bacteria_constants()
    constants.set_simulation_constants()
    constants.set_paths(default=True)

    path = start_run(constants)

    plotting(Path(path))

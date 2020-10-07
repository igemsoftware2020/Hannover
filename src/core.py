#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ********************************************************************************************
# imports

# custom libraries
from src.biofilm import Biofilm
from src.constants import Constants
from src.utils import plot_size, plot_force, plot_velocities, plot_positions, bacteria_as_pandas, \
    prompt_log_at_start, plot_num, dens_map, animate_3d, animate_positions


# ********************************************************************************************
# core function

def blind_run():

    # Set constants for modelling run
    constants = Constants(bac_type="B.Sub.")
    constants.set_bacteria_constants()
    constants.set_simulation_constants()
    constants.set_paths(default=True)

    # Init a new Biofilm with above constants
    biofilm = Biofilm()
    # pass constants to biofilm object
    biofilm.constants = constants

    prompt_log_at_start(biofilm.constants)

    # Save log file for
    info_file_name = constants.get_paths(key="info")
    info_file_path = info_file_name.parent

    biofilm.simulate(save_name=info_file_name)


def plotting(info_file_path):
    data = bacteria_as_pandas(info_file_path)

    plot_num(data, info_file_path, save_fig=True)
    dens_map(data, info_file_path, save_fig=True)
    plot_velocities(data, info_file_path, save_fig=True)
    plot_positions(data, info_file_path, save_fig=True)  # this one messes with data
    plot_force(data, info_file_path, save_fig=True)
    plot_size(data, info_file_path, save_fig=True)
    data = bacteria_as_pandas(info_file_path)
    animate_positions(data, info_file_path, save_fig=True)
    animate_3d(data, info_file_path, save_fig=True)
    
# ********************************************************************************************
# main-method to start the program
# ********************************************************************************************


if __name__ == "__main__":
    blind_run()
    # coreFunction())
    # path = '/home/david/PycharmProjects/biofilm_growth_modeling/output/log_12h59min_2592020/log_12h59min_2592020.json'
    # plotting(path)

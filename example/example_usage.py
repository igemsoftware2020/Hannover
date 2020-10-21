#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ********************************************************************************************
# custom libraries
from BiofilmSimulation.biofilm import Biofilm
from BiofilmSimulation.constants import Constants
from BiofilmSimulation.data_handling import bacteria_as_pandas
from BiofilmSimulation.plotting import histo_length, histo_velocity, histo_force
from BiofilmSimulation.plotting import plot_sizes, plot_force, plot_velocities, plot_positions, \
    animate_positions, plot_num, dens_map
from BiofilmSimulation.utils import prompt_log_at_start


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

    info_file_path = biofilm.constants.get_paths(key="info")
    biofilm.simulate_multiprocessing()
    plotting(info_file_path)


def plotting(info_file_path):
    """ reads in data from info_file_path and plots data """
    data = bacteria_as_pandas(info_file_path)
    histo_length(data, info_file_path, save_fig=True)
    plot_num(data, info_file_path, save_fig=True)
    histo_velocity(data, info_file_path, save_fig=True)
    dens_map(data, info_file_path, save_fig=True)
    plot_velocities(data, info_file_path, save_fig=True)
    histo_force(data, info_file_path, save_fig=True)
    plot_positions(data, info_file_path, save_fig=True)
    plot_force(data, info_file_path, save_fig=True)
    plot_sizes(data, info_file_path, save_fig=True)
    data = bacteria_as_pandas(info_file_path)
    animate_positions(data, info_file_path, save_fig=True)
    # animate_3d(data, info_file_path, save_fig=False)


def default_run():
    strain = input("Select bacteria strain (B.Sub. or E.Coli.) : ")
    num_initial = int(input("Select number of initial bacteria : "))
    duration = int(input("Specify simulation duration in minutes : "))
    time_step = int(input("Set simulation time step in minutes : "))
    constants = Constants(bac_type=strain)
    constants.num_initial_bac = num_initial
    constants.duration = duration
    constants.time_step = time_step
    constants.window_size = (2000, 2000)
    constants.set_bacteria_constants()
    constants.set_simulation_constants()
    constants.set_paths()
    start_run(constants)


# ********************************************************************************************
# main-method to start the program
# ********************************************************************************************


if __name__ == "__main__":
    default_run()

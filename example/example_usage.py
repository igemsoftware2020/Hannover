#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ********************************************************************************************
# custom libraries
from BiofilmSimulation.biofilm import Biofilm
from BiofilmSimulation.constants import Constants
from BiofilmSimulation.data_handling import bacteria_as_pandas, read_in_log, ask_for_log_dir
from BiofilmSimulation.plotting import histo_length, histo_velocity, histo_force
from BiofilmSimulation.plotting import plot_sizes, plot_force, plot_velocities, plot_positions, \
    plot_num, dens_map, animate_positions
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
    biofilm.spawn()
    print(biofilm)


    #print("In Grid", sum(entry is not None for entry in biofilm.bacteria_grid[0]))
    #check_neighbors(biofilm.bacteria_grid, biofilm.coordinates_grid, biofilm.bacteria[0])

    biofilm.simulate_multiprocessing()
    #plotting(info_file_path)


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
    # data = bacteria_as_pandas(info_file_path)
    # animate_positions(data, info_file_path, time_step=time_step, save_fig=True)


def default_run():
    """ 
    Use this function to start a simulation with default values. 
    Enter some simulation paramteres in the terminal, set a reasonable duration and time step
    and automatic plotting results. 
    Keep in mind: Depending on the selected duration, the plot ranges may be not adequate.
    """
    
    # strain = input("Select bacteria strain (B.Sub. or E.Coli.) : ")
    # num_initial = int(input("Select number of initial bacteria : "))
    # duration = int(input("Specify simulation duration in minutes : "))
    # time_step = int(input("Set simulation time step in seconds : "))
    constants = Constants(bac_type="B.Sub.")
    num_initial = 100
    constants.num_initial_bac = num_initial
    #constants.duration = duration
    #constants.time_step = time_step
    constants.set_bacteria_constants()
    constants.set_simulation_constants()
    constants.set_paths()
    start_run(constants)


# ********************************************************************************************
# main-method to start the program
# ********************************************************************************************


if __name__ == "__main__":
    default_run()

    # Use this two commands, if you want to plot data for a specific log
    # path = ask_for_log_dir()
    # plotting(path)


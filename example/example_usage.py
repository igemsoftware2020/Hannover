#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ********************************************************************************************
import argparse
from pathlib import Path

# custom libraries
from BiofilmSimulation.biofilm import Biofilm
from BiofilmSimulation.constants import Constants
from BiofilmSimulation.data_handling import bacteria_as_pandas, read_in_log, ask_for_log_dir
from BiofilmSimulation.plotting import histo_length, histo_velocity, histo_force
from BiofilmSimulation.plotting import plot_sizes, plot_force, plot_velocities, plot_positions, \
    plot_num, dens_map, animate_positions, animate_3d, lennard_jones_force_plot, plot_convex_hull,\
    scatter_last_positions, plot_delauny_triangulation, plot_surface
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
    biofilm.simulate_multiprocessing(center=[500, 500])
    plotting(info_file_path)


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

    #plot_surface(data)
    #plot_convex_hull(data)
    #plot_delauny_triangulation(data)
    #scatter_last_positions(data, info_file_path, save_fig=True)
    #lennard_jones_force_plot(1, 6E-9)
    # Animations
    data = bacteria_as_pandas(info_file_path)
    animate_positions(data, info_file_path, time_step=time_step, save_fig=True)
    animate_3d(data, info_file_path, time_step=time_step, save_fig=False)


def default_run(strain: str = 'E.Coli.', number_initial_bacteria: int = 10, duration: int = 10, time_step: int = 1):
    """ 
    Use this function to start a simulation with default values. 
    Enter some simulation paramteres in the terminal, set a reasonable duration and time step
    and automatic plotting results. 
    Keep in mind: Depending on the selected duration, the plot ranges may be not adequate.
    """
    constants = Constants(bac_type=strain)
    constants.num_initial_bac = number_initial_bacteria
    constants.duration = duration
    constants.time_step = time_step
    constants.set_bacteria_constants()
    constants.set_simulation_constants()
    constants.set_paths()
    start_run(constants)


# ********************************************************************************************
# main-method to start the program
# ********************************************************************************************


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Start Biofilm simulation with custom parameters.')

    parser.add_argument('--strain', type=str, help="Select between 'B.Sub.' or 'E.Coli'.")
    parser.add_argument('-nI', '--numInitial', type=int, help="Specify number of initial Bacteria on the surface.")
    parser.add_argument('-d', '--duration', type=int, help='Duration of simulated time in minutes.')
    parser.add_argument('--step', type=int, help='Specify time steps of the simulation.')
    parser.add_argument('--custom_plot', type=bool,
                        help='If True ist passed, you will be asked for a log file for plotting')

    args = parser.parse_args()

    default_run(strain=args.strain, number_initial_bacteria=args.numInitial,
                duration=args.duration, time_step=args.step)

    # Use this two commands, if you want to plot data for a specific log
    if args.custom_plot:
        path = ask_for_log_dir()
        if path:
            plotting(path)
        else:
            print("No path selected.")


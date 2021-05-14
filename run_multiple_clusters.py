#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ********************************************************************************************
import argparse
import json
from pathlib import Path
# custom libraries
from BiofilmSimulation.constants import Constants
from BiofilmSimulation.biofilm import Biofilm
from BiofilmSimulation.data_handling import read_in_log, bacteria_as_pandas, ask_for_log_dir
from BiofilmSimulation.utils import prompt_log_at_start
from BiofilmSimulation.plotting import plot_sizes, plot_force, plot_velocities, plot_positions, \
    plot_num, dens_map, animate_positions, histo_length, histo_velocity, histo_force


def simulate_n_cluster(constant, centers: [(), (), ()]):
    save_paths = [Path(str(constant.info_path).replace(".json", f"_{n}.json")) for n in range(0, len(centers))]
    for n in range(0, len(centers)):
        biofilm = Biofilm()
        constant.info_path = save_paths[n]
        biofilm.constants = constant
        prompt_log_at_start(biofilm.constants)
        biofilm.center = centers[n]
        biofilm.simulate_multiprocessing()
    return save_paths


def combine_info_files(save_fp: Path, source_fp: [Path]):
    data_unmerged = []
    constants_unmerged = []

    for fp in source_fp:
        data = read_in_log(fp)
        data_unmerged.append(data['BACTERIA'])
        constants_unmerged.append(data['CONSTANTS'])

    # assert all(constants_unmerged[i]['time_step'] == constants_unmerged[0]['time_step']
    #           & constants_unmerged[i]['duration'] == constants_unmerged[0]['duration']
    #           & constants_unmerged[i]['bacteria_strain'] == constants_unmerged[0]['bacteria_strain']
    #           for i in range(0, len(constants_unmerged)))
    merged = {}
    total_count = 0
    for data in data_unmerged:
        count = len(data)
        total_count += count
        if (total_count - count) == 0:
            merged.update(data)
        else:
            for i in range(0, count):
                try:
                    data[f'bacteria_{total_count + i}'] = data.pop(f'bacteria_{i}')
                except KeyError:
                    total_count -= 1
            merged.update(data)

    #assert len(merged) == total_count

    constants_sp = save_fp / 'merged_bacteria_Constants.json'
    bacteria_sp = save_fp / 'merged_bacteria.json'
    with open(constants_sp, 'w+') as fp:
        json.dump(constants_unmerged[0], fp)

    with open(bacteria_sp, 'w+') as fp:
        json.dump({'BACTERIA': merged}, fp)

    return bacteria_as_pandas(bacteria_sp)


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

    constants = Constants(bac_type=args.strain)
    constants.num_initial_bac = args.numInitial
    constants.duration = args.duration
    constants.time_step = args.step
    constants.set_bacteria_constants()
    constants.set_simulation_constants()
    constants.set_paths()

    sources = simulate_n_cluster(constants, centers=[(200, 200), (1000, 200), (200, 1000), (1000, 1000), (400, 400), (600, 750),
                                                     (500, 350), (330, 620), (800, 450)])

    save_fp = (sources[0]).parent
    data = combine_info_files(save_fp=save_fp, source_fp=sources)

    if args.custom_plot:
        path = ask_for_log_dir()
    if path:
        plotting(path)
    else:
        print("No path selected.")


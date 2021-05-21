#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ********************************************************************************************
import argparse
import json
import numpy as np
from pathlib import Path
# custom libraries
from BiofilmSimulation.constants import Constants
from BiofilmSimulation.biofilm import Biofilm
from BiofilmSimulation.data_handling import read_in_log, bacteria_as_pandas, ask_for_log_dir, combine_info_files
from BiofilmSimulation.utils import prompt_log_at_start


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


if __name__ == "__main__":
    # ********************************************************************************************
    # main-method to start the program
    # ********************************************************************************************

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
    x_extend, y_extend, _ = constants.window_size

    x_centers = np.linspace(500, x_extend - 100, 8)
    y_centers = np.linspace(500, y_extend - 100, 8)
    centers = [(xc, yc) for xc in x_centers for yc in y_centers]
    print(centers)
    sources = simulate_n_cluster(constants, centers=centers)
    save_fp = (sources[0]).parent
    data = combine_info_files(save_fp=save_fp, source_fp=sources)

#!/usr/bin/env python
# -*- coding: utf-8 -*-


import random
from pathlib import Path

# ********************************************************************************************
# imports
import numpy as np
import tqdm

# custom libraries
from src.bacteria import Bacterium, get_bacteria_dict
from src.constants import Constants as C
from src.utils import write_log_template, read_in_log, save_dict_as_json, simulation_duration


class Biofilm(object):

    def __init__(self):
        self.bacteria = []
        self.num_bacteria = len(self.bacteria)

    def spawn(self):
        for _ in range(C.NUM_INITIAL_BACTERIA):
            # place bacteria randomly on plate with dimensions C.WINDOW_SIZE[0] um x C.WINDOW_SIZE[1]
            rnd_position = np.asarray([random.uniform(C.WINDOW_SIZE[0] / 2 - 40, C.WINDOW_SIZE[0] / 2 + 40),
                                       random.uniform(C.WINDOW_SIZE[1] / 2 - 40, C.WINDOW_SIZE[1] / 2 + 40), 0])
            # set initial velocity to 0
            velocity = np.asarray([0, 0, 0])
            # random orientation
            rnd_angle = np.asarray([random.uniform(0, 2 * np.pi), random.uniform(0, 2 * np.pi), 0])
            # substrate cell adhesion, in cartesian coordinates
            adhesion_force = np.asarray([0, 0, C.MAX_CELL_SUBSTRATE_ADHESION])
            bac = Bacterium(position=rnd_position, velocity=velocity, angle=rnd_angle, force=adhesion_force,
                            attached_to_surface=True)
            self.bacteria.append(bac)

    def write_to_log(self, log_name):
        info_file_path = C.OUTPUT_PATH / log_name
        if not log_name.is_file():
            # create json template and save
            write_log_template(log_name)

        # read in current log
        data = read_in_log(info_file_path)
        # init empty dictionary. This will overwrite the old entry.
        bacteria_dic = {}
        if sum(map(len, data['BACTERIA'].keys())) == 0:
            # if no entry in BACTERIA, create the first one.
            for bacteria, counter in zip(self.bacteria, range(0, len(self.bacteria))):
                bacteria_dic.update({'bacteria_%s' % str(counter): get_bacteria_dict(bacteria)})

        else:
            # copy already existing one and add to entries
            bacteria_dic = data['BACTERIA'].copy()
            for bacteria, counter in zip(self.bacteria, range(0, len(self.bacteria))):
                bacteria_name = 'bacteria_%s' % str(counter)
                # iterate over all entries in BACTERIA, append next iteration step to key values
                if bacteria_name not in bacteria_dic.keys():
                    # Add bacteria to BACTERIUM keys, because it's not in there
                    bacteria_dic.update({'bacteria_%s' % str(counter): get_bacteria_dict(bacteria)})
                else:
                    for entry in data['BACTERIA'][bacteria_name].keys():
                        # If entry already exists : Append info from next iteration step to corresponding entry
                        for attr in dir(bacteria):
                            if not callable(getattr(bacteria, attr)) and not attr.startswith("__") and entry == attr:
                                attribute = getattr(bacteria, entry)
                                if isinstance(attribute, np.ndarray):
                                    attribute = attribute.tolist()
                                bacteria_dic[bacteria_name][entry].append(attribute)

            # Maybe for checking integrity : len(data['BACTERIA']) - sum(map(len, data['BACTERIA'].keys()))
        data['BACTERIA'] = bacteria_dic
        save_dict_as_json(data, info_file_path)

    @simulation_duration
    def simulate(self, duration_in_min: int, save_name: Path):
        """
        SPAWNS BACTERIA
        Iterates over all Bacteria and updates Parameters
        Because doubling time of bacteria is high with 20 mi, time is measured in minutes.
         All units are SI, therefore conversion in seconds with factor 60 for calculations
        """
        self.spawn()
        print(self)
        print("\nSTARTING MODELLING")
        print(f"SIMULATION TIME INTERVAL {duration_in_min} min in steps of {C.TIME_STEP} s.")
        for _ in tqdm.tqdm(range(0, round(duration_in_min * 60 / C.TIME_STEP))):

            for bacterium in self.bacteria:
                # Grow Bacterium
                bacterium.grow()
                # Split Bacterium
                if bacterium.is_split_ready() and bacterium.living:
                    daughter = bacterium.split()
                    self.bacteria.append(daughter)

                # Forces on bacterium because of drag, adhesion force, gravity
                bacterium.update_acting_force()
                # Add cell- cell interaction force, based on soft-repulsive potential
                for _bacterium in self.bacteria:
                    # check if bacterium is not itself and distance is smaller than 2 times the bacterium length
                    if bacterium != _bacterium and (
                            np.linalg.norm(Biofilm.distance_vector(bacterium, _bacterium)) < 2 * bacterium.length):
                        # add interaction force
                        # bacterium.force += Biofilm.cell_cell_interaction(bacterium, _bacterium)
                        pass

                if bacterium.moving is True:
                    bacterium.update_velocity()
                    bacterium.update_position()
                else:
                    bacterium.moving = random.choices([True, False], weights=[0.5, 0.5])[0]

                if bacterium.living:
                    bacterium.random_cell_death()
                if not bacterium.living:
                    # add increase overall LPS concentration
                    pass

            self.write_to_log(log_name=save_name)

    @staticmethod
    def distance_vector(self: Bacterium, other: Bacterium):
        return self.position - other.position

    @staticmethod
    def cell_cell_interaction(self: Bacterium, other: Bacterium):
        """
        returns force vector of cell- cell interaction.
        Force direction is in direction of the distance vector between the bacteria.
        Force value based on Lennard-Jones Potential / Soft-repulsive potential
        """
        return Biofilm.abs_force_lennard_jones_potential(self, other) * Biofilm.distance_vector(self, other)

    def __repr__(self):
        return f'Biofilm currently consisting of {len(self.bacteria)} bacteria'

    def sort_by_depth(self, axis, _reverse):
        sorted_bacteria = self.bacteria
        # To return a new list, use the sorted() built-in function...
        return sorted(sorted_bacteria, key=lambda x: x.position[axis], reverse=_reverse)

    @staticmethod
    def abs_force_lennard_jones_potential(bacterium1: Bacterium, bacterium2: Bacterium):
        """ return interaction term with bacteria in local environment"""

        def lennard_jones_force(r, epsilon, r_min):
            return - epsilon * (12 * (r_min ** 12 / r ** 13) - 12 * (r_min ** 6 / r ** 7))

        distance_vector = bacterium1.position - bacterium2.position
        distance = np.linalg.norm(distance_vector) * 1E6
        repulsive_force = lennard_jones_force(distance, epsilon=C.MAX_CELL_CELL_ADHESION, r_min=bacterium1.width)
        return repulsive_force

    @staticmethod
    def check_energy_conservation(bacterium1: Bacterium, bacterium2: Bacterium, total_energy_before):
        if bacterium1.total_energy + bacterium2.total_energy != total_energy_before:
            raise ValueError("Energy conversation broken while splitting.")

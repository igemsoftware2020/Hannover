#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
# ********************************************************************************************
# imports
import tqdm

# custom libraries
from src.bacteria import Bacterium, get_bacteria_dict
from src.constants import Constants
from src.utils import write_log_template, read_in_log, save_dict_as_json, simulation_duration


class Biofilm(object):
    """
    The Biofilm class initializes a configuration of bacteria spread on a plane.
    All bacteria of the biofilm are stored as a list and are updated,
    when simulating the growth of the biofilm.
    """

    def __init__(self):
        self.bacteria = []
        self.num_bacteria = len(self.bacteria)
        self.constants = Constants(bac_type="B.Sub.")

    def spawn(self):
        """ spawn an initial number of bacteria.
         Bacteria are randomly distributed on a plane with aspect ratios specified in the Constants class
         """
        num_initial_bacteria = self.constants.get_simulation_constants(key="num_initial")
        window_size = self.constants.get_simulation_constants(key="window_size")
        mean_speed = self.constants.get_bac_constants(key="FREE_MEAN_SPEED")
        for _ in range(num_initial_bacteria):
            # place bacteria randomly on plate with dimensions C.WINDOW_SIZE[0] um x C.WINDOW_SIZE[1]
            rnd_position = np.asarray([np.random.randint(300, 400),
                                       np.random.randint(300, 400),
                                       np.random.normal(3, 0.5)
                                       ])
            # set random initial velocity
            velocity = np.asarray([np.random.normal(mean_speed, mean_speed * 0.01),
                                   np.random.normal(mean_speed, mean_speed * 0.01),
                                   np.random.normal(0, 0.2)
                                   ])
            # random orientation
            rnd_angle = np.asarray([np.random.randint(0, 360),
                                    np.random.normal(0, 360),
                                    np.random.normal(0, 360)
                                    ])
            # substrate cell adhesion, in cartesian coordinates
            adhesion_force = np.asarray([0,
                                         0,
                                         -self.constants.MAX_CELL_SUBSTRATE_ADHESION
                                         ])

            bac = Bacterium(position=rnd_position, velocity=velocity, angle=rnd_angle, force=adhesion_force,
                            attached_to_surface=True, constants=self.constants, strain=self.constants.bac_type)
            self.bacteria.append(bac)

    def write_to_log(self):
        """
        Saves the current parameters of all bacteria in "self.bacteria" as a dictionary in a json file
        with the name log_name. If no json file exits it will create a template. No entries are overwritten,
        instead the parameter lists are updated accordingly
        """
        info_file_path = self.constants.get_paths(key="info")
        if not info_file_path.is_file():
            # create json template and save
            write_log_template(info_file_path, self.constants)

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
    def simulate(self):
        """
        Sets up and runs simulation.
        Iterates over all Bacteria and updates parameters, based on euler- forward scheme.
        All units are SI, expect lengths, which are measured in um
        """
        time_step = self.constants.get_simulation_constants(key="time_step")
        duration = self.constants.get_simulation_constants(key="duration")
        self.spawn()

        print(f"\n ********* STARTING MODELLING  ********* \n "
              f"SIMULATION TIME INTERVAL {duration} min in steps of {time_step} s."
              )
        for _ in tqdm.tqdm(range(0, round(duration * 60 / time_step))):
            cp_bacteria = np.copy(self.bacteria)
            for bacterium in self.bacteria:

                # Grow Bacterium
                bacterium.grow()
                bacterium.update_mass()
                # Split Bacterium
                if bacterium.is_split_ready() and bacterium.living:
                    daughter = bacterium.split()
                    self.bacteria.append(daughter)
                    cp_bacteria = np.append(cp_bacteria, daughter)

                # Forces on bacterium because of drag, adhesion force, gravity
                bacterium.update_acting_force()

                # Add cell- cell interaction force, based on soft-repulsive potential
                for _bacterium in cp_bacteria:
                    # check if bacterium is not itself and distance is smaller than 2 * bacterium length
                    if bacterium != _bacterium \
                            and (np.linalg.norm(Biofilm.distance_vector(bacterium, _bacterium)) < 2 * bacterium.length):
                        # add interaction force
                        force_vector = Biofilm.cell_cell_interaction(bacterium, _bacterium, exact=False)
                        bacterium.force = np.add(bacterium.force, force_vector)
                        _bacterium.force = np.subtract(_bacterium.force, force_vector)
                        # delete bacteria from interaction iteration
                        cp_bacteria = np.delete(cp_bacteria,
                                                np.argwhere((cp_bacteria == _bacterium) & (cp_bacteria == bacterium)))

                bacterium.update_acceleration()
                bacterium.update_velocity()
                bacterium.update_position()

                if bacterium.living is True:
                    bacterium.random_cell_death()
            self.write_to_log()

    @staticmethod
    def distance_vector(self: Bacterium, other: Bacterium):
        """ return distance vector between two bacteria """
        return self.position - other.position

    @staticmethod
    def cell_cell_interaction(self: Bacterium, other: Bacterium, exact=False):
        """
        returns force vector of cell- cell interaction.
        Force direction is in direction of the distance vector between the bacteria.
        Force value based on Lennard-Jones Potential / Soft-repulsive potential
        """
        if exact:
            return Biofilm.abs_force_lennard_jones_potential(bacterium1=self, bacterium2=other, exact=True) \
                   * Biofilm.distance_vector(self, other) / np.linalg.norm(Biofilm.distance_vector(self, other))
        return Biofilm.abs_force_lennard_jones_potential(bacterium1=self, bacterium2=other) * \
               Biofilm.distance_vector(self, other) / np.linalg.norm(Biofilm.distance_vector(self, other))

    def __repr__(self):
        return f'Biofilm consisting of {len(self.bacteria)} bacteria'

    def sort_by_depth(self, axis, _reverse):
        sorted_bacteria = self.bacteria
        # To return a new list, use the sorted() built-in function...
        return sorted(sorted_bacteria, key=lambda x: x.position[axis], reverse=_reverse)

    @staticmethod
    def abs_force_lennard_jones_potential(bacterium1: Bacterium, bacterium2: Bacterium):
        """ return absolute interaction force with one bacteria.
         Interaction force is calculated from the distance and gradient value
         of the lennard- jones potential at this distance
        """

        def lennard_jones_force(r, epsilon, sigma):
            return 48 * epsilon * np.power(sigma, 12) / np.power(r, 13) - 24 * epsilon * np.power(sigma, 6) / np.power(
                r, 7)

        # only calculate for the center of the bacteria
        distance_vector = bacterium1.position - bacterium2.position
        distance = np.linalg.norm(distance_vector) * 1E6
        repulsive_force = lennard_jones_force(distance, epsilon=bacterium1.constants.MAX_CELL_CELL_ADHESION,
                                              sigma=bacterium1.length)
        return repulsive_force

    @staticmethod
    def check_energy_conservation(bacterium1: Bacterium, bacterium2: Bacterium, total_energy_before):
        """ checks if energy is conserved in the splitting process"""
        if bacterium1.total_energy + bacterium2.total_energy != total_energy_before:
            raise ValueError("Energy conversation broken while splitting.")

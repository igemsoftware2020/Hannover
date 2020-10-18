#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ********************************************************************************************
# imports
import numpy as np
import tqdm
from itertools import repeat
from multiprocessing import Pool, cpu_count

# custom libraries
from src.bacteria import Bacterium, get_bacteria_dict
from src.bacteria import distance_vector, bac_bac_interaction_force
from src.constants import Constants
from src.data_handling import write_log_template, read_in_log, save_dict_as_json
from src.utils import simulation_duration


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

    def __repr__(self):
        return f'Biofilm consisting of {len(self.bacteria)} bacteria'

    def spawn(self):
        """
        spawn an initial number of bacteria.
         Bacteria are randomly distributed on a plane with aspect ratios specified in the c class
         """
        num_initial_bacteria = self.constants.get_simulation_constants(key="num_initial")
        window_size = self.constants.get_simulation_constants(key="window_size")
        mean_speed = self.constants.get_bac_constants(key="FREE_MEAN_SPEED")
        for _ in range(num_initial_bacteria - 1):
            # place bacteria randomly on plate with dimensions C.WINDOW_SIZE[0] um x C.WINDOW_SIZE[1]
            rnd_position = np.asarray([np.random.randint(200, 400),
                                       np.random.randint(200, 400),
                                       np.random.normal(3, 0.5)
                                       ])
            # set random initial velocity
            # velocity = np.asarray([np.random.normal(mean_speed, mean_speed * 0.01),
            #                       np.random.normal(mean_speed, mean_speed * 0.01),
            #                       np.random.normal(0, 0.2)
            #                       ])
            # random orientation
            velocity = np.asarray([0, 0, 0])
            rnd_angle = np.asarray([np.random.randint(0, 360),
                                    np.random.randint(0, 360),
                                    np.random.randint(0, 360)
                                    ])
            # substrate cell adhesion, in cartesian coordinates
            adhesion_force = np.asarray([0,
                                         0,
                                         -self.constants.MAX_CELL_SUBSTRATE_ADHESION
                                         ])

            bac = Bacterium(position=rnd_position, velocity=velocity, angle=rnd_angle, force=adhesion_force,
                            attached_to_surface=True, constants=self.constants, strain=self.constants.bac_type)
            self.bacteria.append(bac)

    @simulation_duration
    def simulate_multiprocessing(self):

        time_step = self.constants.get_simulation_constants(key="time_step")
        duration = self.constants.get_simulation_constants(key="duration")
        self.spawn()
        num_threads = cpu_count()
        print(f"\n ********* STARTING MODELLING  USING MULTIPROCESSING ********* \n "
              f"SIMULATION TIME INTERVAL {duration} min in steps of {time_step} s.\n"
              f"Using {num_threads} cores."
              )
        with Pool(processes=num_threads) as pool:
            for _ in tqdm.tqdm(range(0, round(duration * 60 / time_step))):
                try:
                    self.bacteria = pool.map(forces_on_bacterium, self.bacteria)
                    cp_bacteria_list = self.bacteria
                    self.bacteria = pool.starmap(bac_bac_interaction, zip(self.bacteria, repeat(cp_bacteria_list)))
                    self.bacteria = pool.map(update_movement, self.bacteria)

                    for mother in self.bacteria:
                        if mother.is_split_ready() and mother.living:
                            daughter = mother.split()
                            self.bacteria.append(daughter)

                    self.bacteria = pool.map(grow_bacterium, self.bacteria)
                    self.write_to_log()

                except KeyboardInterrupt:
                    self.write_to_log()
                    return self.constants.get_paths(key="info")

            self.write_to_log()
            return self.constants.get_paths(key="info")

    def write_to_log(self):
        """
        Saves the current parameters of all bacteria in "self.bacteria" as a dictionary in a json file
        with the name log_name. If no json file exits it will create a template. No entries are overwritten,
        instead the parameter lists are updated accordingly
        The dictionary is build like this and stored as a json
        {
            {BACTERIA: {bacteria_0:
                            {
                                position: [[] , ... ,[]],
                                velocity: [[] , ... ,[]],
                                ...}
                            .
                            .
                            .
                        bacteria_n:
                                {...}
                        },
            CONSTANTS: { ... }
        }
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
        save_dict_as_json(data['CONSTANTS'], str(info_file_path).replace(".json", "_Constants.json"))
        save_dict_as_json(data, info_file_path)


# These functions are for the multithreading simulation
def grow_bacterium(bacterium: Bacterium):
    # Grow Bacterium
    bacterium.grow()
    bacterium.update_mass()
    return bacterium


def forces_on_bacterium(bacterium: Bacterium):
    bacterium.update_acting_force()
    return bacterium


def update_movement(bacterium: Bacterium):
    bacterium.update_acceleration()
    if np.linalg.norm(bacterium.acceleration) > 0.1:
        bacterium.acceleration = bacterium.acceleration / np.linalg.norm(bacterium.acceleration) * np.random.normal(0,
                                                                                                                    scale=0.002)
    bacterium.update_velocity()
    if not np.linalg.norm(bacterium.velocity > 14):
        bacterium.velocity = bacterium.velocity / 2
    bacterium.update_position()

    if bacterium.living is True:
        bacterium.random_cell_death()
        bacterium.detach()

    return bacterium


def bac_bac_interaction(bacterium: Bacterium, bac_list):
    for other_bacterium in bac_list:
        # this is the bacteria i want to update the interaction force on
        if bacterium != other_bacterium \
                and (
                np.linalg.norm(distance_vector(bacterium, other_bacterium)) < 2 * bacterium.length):
            # add interaction force
            force_vector = bac_bac_interaction_force(bacterium, other_bacterium)
            bacterium.force = np.add(bacterium.force, force_vector)
            bacterium.total_force = np.linalg.norm(bacterium.force)
    return bacterium

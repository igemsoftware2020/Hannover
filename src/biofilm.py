#!/usr/bin/env python
# -*- coding: utf-8 -*-


import json
import random
from typing import Dict
from pathlib import Path
import matplotlib.pyplot as plt

# ********************************************************************************************
# imports
import numpy as np
import tqdm
import glob
import pandas as pd

# custom libraries
from src.bacteria import Bacterium
from src.constants import Constants as C


class Biofilm(object):

    def __init__(self):
        self.bacteria = []

    def spawn(self):
        for _ in range(C.START_NUMBER_BACTERIA):
            pos = np.asarray(
                [(random.random() - 0.5) * C.WINDOW_SIZE[0], (random.random() - 0.5) * C.WINDOW_SIZE[1], 0])
            velocity = np.asarray([random.gauss(C.BSUB_MEAN_SPEED, 10), random.gauss(C.BSUB_MEAN_SPEED, 10), 0])
            bac = Bacterium(position=pos, velocity=velocity)
            self.bacteria.append(bac)

    @staticmethod
    def check_energy_conservation(bacterium1: Bacterium, bacterium2: Bacterium, total_energy_before):
        if bacterium1.total_energy + bacterium2.total_energy != total_energy_before:
            raise ValueError("Energy conversation broken while splitting.")

    def evolve(self):
        # Iterate over time steps
        for bacterium in self.bacteria:
            # Iterate over bacteria in self.bacteria
            bacterium.grow()

            if not bacterium.living and bacterium.get_volume() < 3300:
                self.bacteria.remove(bacterium)
            # Manage repulsion
            # (The most expensive task)
            for _bacterium in self.bacteria:
                if not np.array_equal(bacterium.position, _bacterium.position):
                    [_bacterium, bacterium] = Biofilm.interaction(bacterium, _bacterium)

            bacterium.move()

            # Manage Bacterial splitting
            # Add a little bit of random until -----it looks good and real-----
            if bacterium.is_split_ready() and bacterium.living:
                energy_before = bacterium.total_energy
                daughter = bacterium.split()
                # self.check_energy_conservation(bacterium, daughter, energy_before)
                self.bacteria.append(daughter)

            # Manage Bacterial Motion-Mode
            # Enable/Disable motion mode
            if bacterium.living:
                if random.random() > 1.0 - C.MOTION_ACTIVATION_PROBABILITY:
                    bacterium.moving = True
                if random.random() > 1.0 - C.MOTION_DEACTIVATION_PROBABILITY:
                    bacterium.moving = False
            else:
                bacterium.moving = False

    @staticmethod
    def write_log_template(info_file_path):
        constants = C()
        with open(info_file_path, 'w+') as json_file:
            data = {'BACTERIA': {}, 'CONSTANTS': {}}
            members = [attr for attr in dir(constants) if
                       not callable(getattr(constants, attr)) and not attr.startswith("__")]
            constants_dic = {}
            for constant in members:
                constants_dic.update({constant: str(getattr(constants, constant))})
            data['CONSTANTS'].update(constants_dic)
            json.dump(data, json_file, indent=4)

    @staticmethod
    def read_in_log(info_file_path) -> Dict:
        with open(info_file_path, "r") as json_file:
            data = json.load(json_file)
        return data

    @staticmethod
    def bacteria_as_pandas(info_file_path) -> pd.DataFrame:
        data = Biofilm.read_in_log(info_file_path)
        return pd.DataFrame(data['BACTERIA']).transpose()

    @staticmethod
    def constants_as_pandas(info_file_path):
        data = Biofilm.read_in_log(info_file_path)
        return pd.DataFrame(data['CONSTANTS']).transpose()

    @staticmethod
    def save_dict_as_json(data: Dict, info_file_path: Path):
        with open(str(info_file_path), 'w') as json_file:
            json.dump(data, json_file)


    def write_to_log(self, log_name):
        info_file_path = C.OUTPUT_PATH / log_name

        def bacteria_dict(bacterium: Bacterium, number: int) -> Dict:
            """ Help function, which returns the dict entry of a bacteria """
            return {'bacteria_%s' % str(number): dict(
                                                position=[bacterium.position.tolist()],
                                                velocity=[bacterium.velocity.tolist()],
                                                angle=[bacterium.angle],
                                                total_force=[bacterium.total_force],
                                                total_energy=[bacterium.total_energy],
                                                living=[bacterium.living],
                                                moving=[bacterium.moving],
                                                length=[bacterium.length],
                                                width=[bacterium.width])
                    }

        if not info_file_path.is_file():
            # create json template and save
            self.write_log_template(info_file_path)

        # read in current log
        data = self.read_in_log(info_file_path)
        # init empty dictionary. This will overwrite the old entry.
        bacteria_dic = {}
        if sum(map(len, data['BACTERIA'].keys())) == 0:
            # if no entry in BACTERIA, create the first one.
            for bacteria, counter in zip(self.bacteria, range(0, len(self.bacteria))):
                bacteria_dic.update(bacteria_dict(bacteria, counter))

        else:
            # copy already existing one and add to entries
            bacteria_dic = data['BACTERIA'].copy()
            for bacteria, counter in zip(self.bacteria, range(0, len(self.bacteria))):
                bacteria_name = 'bacteria_%s' % str(counter)
                # iterate over all entries in BACTERIA, append next iteration step to key values
                if bacteria_name not in bacteria_dic.keys():
                    # Add bacteria to BACTERIUM keys, because it's not in there
                    bacteria_dic.update(bacteria_dict(bacteria, counter))
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
        self.save_dict_as_json(data, info_file_path)

    @staticmethod
    def interaction(bacterium1: Bacterium, bacterium2: Bacterium):
        """ return interaction term with bacteria in local environment"""

        # TODO: update velocities accordingly
        def lennard_jones_force(r, epsilon, r_min):
            return - epsilon * (12 * (r_min ** 12 / r ** 13) - 12 * (r_min ** 6 / r ** 7))

        distance_vector = bacterium1.position - bacterium2.position
        distance = np.sqrt(np.dot(distance_vector, distance_vector))

        if distance < bacterium1.length * 4 and distance != 0:
            # Values for epsilon and r_min kinda random
            repulsive_force = lennard_jones_force(distance, epsilon=0.1, r_min=bacterium1.length)
            acceleration = repulsive_force / C.BSUB_MASS

            def new_position(a, v, t, r): return 1 / 2 * a * t ** 2 + v * t + r

            bacterium1.position = new_position(acceleration, bacterium1.velocity, C.TIME_STEP, bacterium1.position)
            bacterium2.position = new_position(acceleration, bacterium2.velocity, C.TIME_STEP, bacterium2.position)

        return bacterium1, bacterium2

    def __repr__(self):
        return f'Biofilm currently consisting of {len(self.bacteria)} bacteria'

    def sort_by_depth(self, axis, _reverse):
        sorted_bacteria = self.bacteria
        # To return a new list, use the sorted() built-in function...
        return sorted(sorted_bacteria, key=lambda x: x.position[axis], reverse=_reverse)

    @staticmethod
    def plot_velocities(data: pd.DataFrame):
        plot_data = []
        for index, bacteria in data.iterrows():
            velocities = pd.DataFrame(bacteria).loc['velocity',:]
            velocities = velocities.transform(lambda x: sorted(x, key=pd.isnull, reverse=True))
            for vectors in velocities:
                vectors = (pd.Series(vectors).apply(np.array)).apply(np.linalg.norm)
                plot_data.append(vectors)

        for data in plot_data:
            plt.plot(data)
        # plt.plot(, label='mean')  # plot means for each iteration
        plt.title(' VELOCITIES')
        plt.legend()
        plt.show()


    @staticmethod
    def plot_xy_trajectories(data):
        fig, ax = plt.subplots()
        x_data, y_data = [], []
        for bacteria in data['BACTERIA'].keys():
            x, y = data['BACTERIA'][bacteria]['position'][0], data['BACTERIA'][bacteria]['position'][1]
            x_data.append(x)
            y_data.append(y)

        for i in range(0, len(x_data)):
            ax.plot(x_data[i], y_data[i], '.')
            
        plt.show()
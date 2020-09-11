#!/usr/bin/env python
# -*- coding: utf-8 -*-


import random

# ********************************************************************************************
# imports
import numpy as np
import tqdm

# custom libraries
from src.bacteria import Bacterium
from src.constants import Constants as C
from src.utils import get_info_file_path
from src.utils import write_log_template, read_in_log, save_dict_as_json


class Biofilm(object):

    def __init__(self):
        self.bacteria = []
        self.num_bacteria = len(self.bacteria)

    def spawn(self):
        for _ in range(C.START_NUMBER_BACTERIA):
            # place bacteria randomly on plate with dimensions C.WINDOW_SIZE[0] um x C.WINDOW_SIZE[1]
            rnd_position = np.asarray([random.randrange(C.WINDOW_SIZE[0] / 2 - 40, C.WINDOW_SIZE[0] / 2 + 40),
                                       random.randrange(C.WINDOW_SIZE[1] / 2 - 40, C.WINDOW_SIZE[1] / 2 + 40), 0])
            # set initial velocity to 0
            velocity = np.asarray([0, 0, 0])
            # random orientation
            rnd_angle = np.ndarray([random.randrange(0, 2 * np.pi), random.randrange(0, 2 * np.pi)])
            # substrate cell adhesion, in cartesian coordinates
            adhesion_force = np.asarray([0, 0, C.MAX_CELL_SUBSTRATE_ADHESION])
            bac = Bacterium(position=rnd_position, velocity=velocity, angle=rnd_angle, force=adhesion_force)
            self.bacteria.append(bac)

    def write_to_log(self, log_name):
        info_file_path = C.OUTPUT_PATH / log_name
        if not info_file_path.is_file():
            # create json template and save
            write_log_template(info_file_path)

        # read in current log
        data = read_in_log(info_file_path)
        # init empty dictionary. This will overwrite the old entry.
        bacteria_dic = {}
        if sum(map(len, data['BACTERIA'].keys())) == 0:
            # if no entry in BACTERIA, create the first one.
            for bacteria, counter in zip(self.bacteria, range(0, len(self.bacteria))):
                bacteria_dic = bacteria.get_bacteria_dict()
                bacteria_dic.update({'bacteria_%s' % str(counter): bacteria_dic})

        else:
            # copy already existing one and add to entries
            bacteria_dic = data['BACTERIA'].copy()
            for bacteria, counter in zip(self.bacteria, range(0, len(self.bacteria))):
                bacteria_name = 'bacteria_%s' % str(counter)
                # iterate over all entries in BACTERIA, append next iteration step to key values
                if bacteria_name not in bacteria_dic.keys():
                    # Add bacteria to BACTERIUM keys, because it's not in there
                    bacteria_dic = bacteria.get_bacteria_dict()
                    bacteria_dic.update({'bacteria_%s' % str(counter): bacteria_dic})
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

    def simulate(self, duration_in_min: int):
        """
        SPAWNS BACTERIA
        Iterates over all Bacteria and updates Parameters
        Because doubling time of bacteria is high with 20 mi, time is measured in minutes.
         All units are SI, therefore conversion in seconds with factor 60 for calculations
        """
        info_file_name = get_info_file_path()
        self.spawn()
        print(self)
        print("\nSTARTING MODELLING")
        print(f"SIMULATE TIME INTERVAL {duration_in_min} min in steps of {C.TIME_STEP} s.")
        for _ in tqdm.tqdm(range(0, duration_in_min * 60 / C.TIME_STEP)):
            for bacterium in self.bacteria:
                # Grow Bacterium
                bacterium.grow()
                # Split Bacterium
                if bacterium.is_split_ready() and bacterium.living:
                    daughter = bacterium.split()
                    self.bacteria.append(daughter)
                # Forces on bacterium because of drag, adhesion force, gravity
                bacterium.update_acting_force()
                for _bacterium in self.bacteria:
                    if _bacterium != bacterium and not np.array_equal(bacterium.position, _bacterium.position) and \
                            (np.linalg.norm(Biofilm.distance_vector(bacterium, _bacterium)) < 2 * bacterium.length):
                        bacterium.force += Biofilm.cell_cell_interaction(bacterium, _bacterium)

                bacterium.update_velocity()
                bacterium.move()

                bacterium.random_cell_death()

            self.write_to_log(log_name=info_file_name)

    @staticmethod
    def distance_vector(self: Bacterium, other: Bacterium):
        return self.position - other.position

    @staticmethod
    def cell_cell_interaction(self: Bacterium, other: Bacterium):
        return Biofilm.abs_force_lennard_jones_potential(other, other) * Biofilm.distance_vector(self, other)

    def __repr__(self):
        return f'Biofilm currently consisting of {len(self.bacteria)} bacteria'

    def sort_by_depth(self, axis, _reverse):
        sorted_bacteria = self.bacteria
        # To return a new list, use the sorted() built-in function...
        return sorted(sorted_bacteria, key=lambda x: x.position[axis], reverse=_reverse)

    @staticmethod
    def abs_force_lennard_jones_potential(bacterium1: Bacterium, bacterium2: Bacterium):
        """ return interaction term with bacteria in local environment"""

        # TODO: update velocities accordingly
        def lennard_jones_force(r, epsilon, r_min):
            return - epsilon * (12 * (r_min ** 12 / r ** 13) - 12 * (r_min ** 6 / r ** 7))

        distance_vector = bacterium1.position - bacterium2.position
        distance = np.linalg.norm(distance_vector)
        repulsive_force = lennard_jones_force(distance, epsilon=C.MAX_CELL_CELL_ADHESION, r_min=bacterium1.width)
        return repulsive_force

    @staticmethod
    def interaction(bacterium1, bacterium2):
        """ return interaction term with bacteria in local environment"""
        lenPos = len(bacterium1.get_position())

        dx = bacterium2.position[0] - bacterium1.position[0]
        dy = bacterium2.position[1] - bacterium1.position[1]
        dz = bacterium2.position[2] - bacterium1.position[2]
        dr = dx * dx + dy * dy + dz * dz
        interactionfactor = 0.5
        # If the bacterium is "far away"
        # Ignore all circles of this Bacteria to improve speed
        # Do the following operation with all other Bacteria
        far_away_factor = 8
        if (dr ** 0.5 < (bacterium1.width) * 1.05 * 5 * far_away_factor):
            positions = bacterium1.get_position()
            for index in range(lenPos - 1):
                position = positions[index]

                dx = bacterium2.position[0] - position[0]  # self.pos[0]
                dy = bacterium2.position[1] - position[1]  # self.pos[1]
                dz = bacterium2.position[2] - position[2]  # self.pos[2]
                dr = dx * dx + dy * dy + dz * dz
                interactionfactor = 0.5
                # If the bacterium is "far away"
                # Ignore all circles of this one to improve speed
                far_away_factor = 4
                if (dr ** 0.5 < (bacterium1.width) * 1.05 * 5 * far_away_factor):

                    # Each not only with the center of _Bacterium, instead also every sphere of this
                    _positions = bacterium2.get_position()
                    _lenPos = len(_positions)
                    for _index in range(_lenPos - 1):
                        _position = _positions[_index]

                        dx = _position[0] - position[0]  # self.pos[0]
                        dy = _position[1] - position[1]  # self.pos[1]
                        dz = _position[2] - position[2]  # self.pos[2]
                        dr = dx * dx + dy * dy + dz * dz
                        interactionfactor = 0.25
                        # if True:

                        if (dx != 0) or (dy != 0) or (dz != 0):

                            repulsion_x = 0
                            repulsion_y = 0
                            repulsion_z = 0

                            # if(dr**0.5<(Bacterium.getVolume())**(1/3)*1.25*1.65):
                            #    repulsion_x = -dx*150    /dr**(1.5)
                            #    repulsion_y = -dy*150    /dr**(1.5)
                            if (dr ** 0.5 < (bacterium1.width) * 1.05 * 5):
                                repulsion_x = -dx * 40 / dr  # (1.5)
                                repulsion_y = -dy * 40 / dr  # (1.5)
                                repulsion_z = -dz * 40 / dr  # (1.5)

                            # New repulsion-function design
                            # if(dr**0.5<(Bacterium.getVolume())**(1/3)*1.25*1.7):
                            #    repulsion_x = -dx*40    /dr**(0.5)*repulsion_function(Bacterium.getVolume(), dr)
                            #    repulsion_y = -dy*40    /dr**(0.5)*repulsion_function(Bacterium.getVolume(), dr)

                            #
                            bacterium1.velocity[0] = bacterium1.velocity[0] + int(
                                repulsion_x / lenPos ** (1 / 2) * interactionfactor)
                            bacterium1.velocity[1] = bacterium1.velocity[1] + int(
                                repulsion_y / lenPos ** (1 / 2) * interactionfactor)
                            bacterium1.velocity[2] = bacterium1.velocity[2] + int(
                                repulsion_z / lenPos ** (1 / 2) * interactionfactor)
                            # add torque
                            t_radius = (index - lenPos * 0.5)
                            bacterium1.velocity_angular[0] = bacterium1.velocity_angular[0] + t_radius * np.cos(
                                bacterium1.angle[0]) * np.cos(
                                bacterium1.angle[1]) * repulsion_x / lenPos * 0.05 * interactionfactor
                            bacterium1.velocity_angular[0] = bacterium1.velocity_angular[0] - t_radius * np.sin(
                                bacterium1.angle[0]) * np.cos(
                                bacterium1.angle[1]) * repulsion_y / lenPos * 0.05 * interactionfactor
                            # Torque on second angle
                            bacterium1.velocity_angular[1] = bacterium1.velocity_angular[1] + t_radius * np.cos(bacterium1.angle[1]) * (
                                        repulsion_x ** 2 + repulsion_y ** 2) ** (
                                                                   1 / 2) / lenPos * 0.0125 * interactionfactor
                            bacterium1.velocity_angular[1] = bacterium1.velocity_angular[1] + t_radius * np.sin(
                                bacterium1.angle[1]) * repulsion_z / lenPos * 0.05 * interactionfactor
                            bacterium1.total_force = bacterium1.total_force + np.linalg.norm(
                                repulsion_x / lenPos * interactionfactor)
                            bacterium1.total_force = bacterium1.total_force + np.linalg.norm(
                                repulsion_y / lenPos * interactionfactor)

                            # Actio-Reactio
                            bacterium2.velocity[0] = bacterium2.velocity[0] - (
                                        repulsion_x / lenPos ** (1 / 2) * interactionfactor)
                            bacterium2.velocity[1] = bacterium2.velocity[1] - (
                                        repulsion_y / lenPos ** (1 / 2) * interactionfactor)
                            bacterium2.velocity[2] = bacterium2.velocity[2] - (
                                        repulsion_z / lenPos ** (1 / 2) * interactionfactor)
                            bacterium2.total_force = bacterium2.total_force + np.linalg.norm(
                                repulsion_y / lenPos * interactionfactor)
                            bacterium2.total_force = bacterium2.total_force + np.linalg.norm(
                                repulsion_x / lenPos * interactionfactor)

        return bacterium2, bacterium1

    @staticmethod
    def check_energy_conservation(bacterium1: Bacterium, bacterium2: Bacterium, total_energy_before):
        if bacterium1.total_energy + bacterium2.total_energy != total_energy_before:
            raise ValueError("Energy conversation broken while splitting.")

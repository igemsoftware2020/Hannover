#!/usr/bin/env python
# -*- coding: utf-8 -*-


# ********************************************************************************************
# imports
import random
import math
import numpy as np
from src.constants import Constants as C


# ********************************************************************************************
# main class <<<bacterium>>>


class Bacterium:

    def __init__(self, position: np.ndarray = None,
                 width: float = C.BSUB_WIDTH,
                 length: float = C.BSUB_LENGTH,
                 velocity: np.ndarray = None,
                 angle=None, total_force: float = 0,
                 living: bool = True, moving: bool = True):
        """
        initialize a instance of the Bacteria class
        :param position: position of bacteria center [x_pox, y_pos, z_pos]
        :param width: width of ellipse in meter
        :param length:  length of ellipse in meter, default value 2 µm for B. subtilis
        :param velocity: velocity of bacteria [v_x, v_y, v_z] in m/s
        :param angle: angle of bacteria  measured to x axis in radian

        """

        if angle is None:
            angle = [random.randint(0,360), random.randint(0,360)]
        else:
            self.angle = angle
        if position is None:
            self.position = [C.WINDOW_SIZE[0] / 2, C.WINDOW_SIZE / 2, 0]
        else:
            self.position = position
        if velocity is None:
            self.velocity = np.asarray([0, 0, 0], dtype=np.int64)
        else:
            self.velocity = np.asarray([velocity[0], velocity[1], velocity[2]], dtype=np.int64)

        self.width = width
        self.length = length
        self.angle = angle  # angle is orientation of cell measured to x- axes in degree
        self.living = living
        self.moving = moving

        self.velocity_angular = [0, 0]
        self.total_force = total_force
        self.rotational_energy = self.get_rotational_energy()
        self.translation_energy = self.get_translation_energy()
        self.total_energy = self.translation_energy + self.rotational_energy

    def get_rotational_energy(self):
        moment_of_inertia = C.BSUB_MASS / 6 * (3 * self.width ** 2 + self.length) + C.BSUB_MASS / 2 * self.width ** 2
        self.rotational_energy = 1 / 2 * moment_of_inertia * np.dot(self.velocity_angular, self.velocity_angular)
        return self.rotational_energy

    def get_translation_energy(self) -> float:
        return 1 / 2 * C.BSUB_MASS * np.dot(self.velocity, self.velocity)

    def get_volume(self):
        """ gives out the cubic volume equivalent """
        return self.width * self.width * self.length

    def grow(self):
        """
        grow bacteria for 1 second with speed growth_rate
        """
        # Make the bacteria grow
        # using a constant growth rate

        if self.living and self.moving:
            growth_suppressor = C.BSUB_GROWTH_FACTOR * C.gr_pr_i / (C.gr_pr_i + self.total_force * 0.5) \
                                - C.gr_factor_inv
            volume = self.get_volume()
            growth_factor = (volume + volume ** (1 / 3) * 5.0 * growth_suppressor * growth_suppressor) / volume
            # self.width  = self.width *growth_factor
            self.length = self.length * (1 + growth_factor)
        else:
            # self.width  = self.width *gr_d_factor
            # size decreasing if cell is dead
            self.length = self.length * (1 - 0.05)

        # self.random_cell_death()

    def random_cell_death(self):
        # Programmed cell death
        if random.random() > 1.0 - C.BSUB_MORTALITY_RATE:
            self.living = False

    def split(self):
        """
        split bacteria  create new daughter bacteria with new values and update the mother bacterium
        :return: daughter bacteria
        """
        # Calculate new position and angle
        # Advanced new position: random angular component,
        #                       radial component sum of two radii*0.8
        def update_angle(angle):
            return np.array(angle) + (0.5 - random.random()) * np.pi * 0.5

        def get_daughter_position(position, split_distance, angle):
            offset = (split_distance * math.sin(angle[0]) * math.cos(angle[1]),
                      split_distance * math.cos(angle[0]) * math.cos(angle[1]),
                      split_distance * math.sin(angle[1]))
            position = position + np.asarray(offset) * volume_ratio
            return position

        volume_ratio = 0.4 + 0.2 * random.random()
        # Create daughter bacterium from self

        # Update parameters of daughter and mother bacterium
        daughter_bac_length = volume_ratio * self.length
        daughter_bac_angle = update_angle(self.angle)
        daughter_bac_position = get_daughter_position(position=self.position, split_distance=self.length * 0.2,
                                                      angle=daughter_bac_angle)
        daughter_bac = Bacterium(daughter_bac_position, self.width, daughter_bac_length,
                                 self.velocity, self.angle, 0, True, False)

        # update mother cell
        self.length = (1 - volume_ratio) * self.length
        return daughter_bac

    def get_position(self) -> np.ndarray:
        positions = []
        dx_length = self.length * math.sin(self.angle[0]) * math.cos(self.angle[1])
        dy_length = self.length * math.cos(self.angle[0]) * math.cos(self.angle[1])
        dz_length = self.length * math.sin(self.angle[1])
        positions.append((self.position[0] - int(0.1 * dx_length), self.position[1] - int(0.1 * dy_length),
                          self.position[2] - int(0.1 * dz_length)))
        for index in range(-int(self.length * 0.01), int(self.length * 0.01) + 1):
            positions.append((self.position[0] + int(10 * index * math.sin(self.angle[0]) * math.cos(self.angle[1])),
                              self.position[1] + int(10 * index * math.cos(self.angle[0]) * math.cos(self.angle[1])),
                              self.position[2] + int(10 * index * math.sin(self.angle[1]))))
        positions.append((self.position[0] + int(0.1 * dx_length), self.position[1] + int(0.1 * dy_length),
                          self.position[2] + int(0.1 * dz_length)))
        return np.asarray(positions)

    def move(self):
        """
        Move Bacteria for 1 time unit and add random movement and rotation
        :return:
        """
        # Set velocities depending on bacteria state
        if self.moving and self.living:
            self.velocity[0] = (np.linalg.norm(self.velocity) * math.cos(self.angle[0]) * math.cos(self.angle[1])) * (
                    1 + C.TIME_STEP)  # calculates x-vel before turning again
            self.velocity[1] = (np.linalg.norm(self.velocity) * math.sin(self.angle[0]) * math.cos(self.angle[1])) * (
                    1 + C.TIME_STEP)  # calculates x-vel before turning again

            # Davids Code for turning goes here
            self.velocity_angular[0] = self.velocity_angular[0] + (0.5 - random.random()) * 0.1 * C.TIME_STEP
            self.velocity_angular[1] = self.velocity_angular[1] + (0.5 - random.random()) * 0.001 * C.TIME_STEP
            # slight z-brownian random drift
            self.velocity[2] = self.velocity[2] + (0.5 - random.random()) * 0.1 * C.TIME_STEP
            # And gravity
            self.velocity[2] = self.velocity[2] - 9.81 * 0.5 * C.TIME_STEP

        if not self.moving and self.living:
            self.velocity_angular[0] = self.velocity_angular[0] + (0.5 - random.random()) * 0.01 * C.TIME_STEP
            self.velocity_angular[1] = self.velocity_angular[1] + (0.5 - random.random()) * 0.001 * C.TIME_STEP

        if self.at_boundary() == "X":
            self.velocity[0] = - self.velocity[0]
        if self.at_boundary() == "Y":
            self.velocity[1] = - self.velocity[1]

        self.position = self.position + self.velocity

    def at_boundary(self):
        x, y, z = self.position
        if x + self.length >= C.WINDOW_SIZE[0] or x - self.length <= 0:
            # elastic scattering at boundary
            return "X"
        elif y + self.length >= C.WINDOW_SIZE[1] or y - self.length <= 0:
            return "Y"
        return False

    def is_split_ready(self):
        def gaussian_distribution(x, mu=C.BSUB_CRITICAL_LENGTH, sigma=C.BSUB_CRITICAL_LENGTH * 0.12):
            return 1 / np.sqrt(2 * np.pi * sigma ** 2) * np.exp(- (x - mu) ** 2 / (2 * sigma ** 2))

        splitting_lengths = np.random.normal(C.BSUB_CRITICAL_LENGTH, C.BSUB_CRITICAL_LENGTH * 0.12)
        if splitting_lengths <= self.length:
            # probability_to_split = gaussian_distribution(self.length)
            # return np.random.choice([True, False], p=[probability_to_split, 1 - probability_to_split])
            return True
        else:
            return False

    def interaction(self, _bacterium):
        """ return interaction term with bacteria in local environment"""
        len_pos = len(self.get_position())

        dx = _bacterium.pos[0] - self.position[0]
        dy = _bacterium.pos[1] - self.position[1]
        dz = _bacterium.pos[2] - self.position[2]
        dr = dx * dx + dy * dy + dz * dz
        # If the bacterium is "far away"
        # Ignore all circles of this Bacteria to improve speed
        # Do the following operation with all other Bacteria
        far_away_factor = 8
        if dr ** 0.5 < self.width * 1.05 * 5 * far_away_factor:
            positions = self.get_position()
            for index in range(len_pos - 1):
                position = positions[index]

                dx = _bacterium.pos[0] - position[0]  # self.pos[0]
                dy = _bacterium.pos[1] - position[1]  # self.pos[1]
                dz = _bacterium.pos[2] - position[2]  # self.pos[2]
                dr = dx * dx + dy * dy + dz * dz
                # If the bacterium is "far away"
                # Ignore all circles of this one to improve speed
                far_away_factor = 4
                if dr ** 0.5 < self.width * 1.05 * 5 * far_away_factor:

                    # Each not only with the center of _Bacterium, instead also every sphere of this
                    _positions = _bacterium.get_position()
                    _lenPos = len(_positions)
                    for _index in range(_lenPos - 1):
                        _position = _positions[_index]

                        dx = _position[0] - position[0]  # self.pos[0]
                        dy = _position[1] - position[1]  # self.pos[1]
                        dz = _position[2] - position[2]  # self.pos[2]
                        dr = dx * dx + dy * dy + dz * dz
                        interaction_factor = 0.25
                        # if True:

                        if (dx != 0) or (dy != 0) or (dz != 0):

                            repulsion_x = 0
                            repulsion_y = 0
                            repulsion_z = 0

                            # if(dr**0.5<(Bacterium.getVolume())**(1/3)*1.25*1.65):
                            #    repulsion_x = -dx*150    /dr**(1.5)
                            #    repulsion_y = -dy*150    /dr**(1.5)
                            if dr ** 0.5 < self.width * 1.05 * 5:
                                repulsion_x = -dx * 40 / dr  # (1.5)
                                repulsion_y = -dy * 40 / dr  # (1.5)
                                repulsion_z = -dz * 40 / dr  # (1.5)

                            # New repulsion-function design
                            # if(dr**0.5<(Bacterium.getVolume())**(1/3)*1.25*1.7):
                            #    repulsion_x = -dx*40    /dr**(0.5)*repulsion_function(Bacterium.getVolume(), dr)
                            #    repulsion_y = -dy*40    /dr**(0.5)*repulsion_function(Bacterium.getVolume(), dr)

                            #
                            self.velocity[0] = self.velocity[0] + int(
                                repulsion_x / len_pos ** (1 / 2) * interaction_factor)
                            self.velocity[1] = self.velocity[1] + int(
                                repulsion_y / len_pos ** (1 / 2) * interaction_factor)
                            self.velocity[2] = self.velocity[2] + int(
                                repulsion_z / len_pos ** (1 / 2) * interaction_factor)
                            # add torque
                            t_radius = (index - len_pos * 0.5)
                            self.velocity_angular[0] = self.velocity_angular[0] + t_radius * math.cos(
                                self.angle[0]) * math.cos(
                                self.angle[1]) * repulsion_x / len_pos * 0.05 * interaction_factor
                            self.velocity_angular[0] = self.velocity_angular[0] - t_radius * math.sin(
                                self.angle[0]) * math.cos(
                                self.angle[1]) * repulsion_y / len_pos * 0.05 * interaction_factor
                            # Torque on second angle
                            self.velocity_angular[1] = self.velocity_angular[1] + t_radius * math.cos(self.angle[1]) * (
                                        repulsion_x ** 2 + repulsion_y ** 2) ** (
                                                                   1 / 2) / len_pos * 0.0125 * interaction_factor
                            self.velocity_angular[1] = self.velocity_angular[1] + t_radius * math.sin(
                                self.angle[1]) * repulsion_z / len_pos * 0.05 * interaction_factor
                            self.total_force = self.total_force + np.abs(
                                repulsion_x / len_pos * interaction_factor)
                            self.total_force = self.total_force + np.abs(
                                repulsion_y / len_pos * interaction_factor)

                            # Actio-Reactio
                            _bacterium.velocity[0] = _bacterium.velocity[0] - (
                                        repulsion_x / len_pos ** (1 / 2) * interaction_factor)
                            _bacterium.velocity[1] = _bacterium.velocity[1] - (
                                        repulsion_y / len_pos ** (1 / 2) * interaction_factor)
                            _bacterium.velocity[2] = _bacterium.velocity[2] - (
                                        repulsion_z / len_pos ** (1 / 2) * interaction_factor)
                            _bacterium.total_force = _bacterium.total_force + np.abs(
                                repulsion_y / len_pos * interaction_factor)
                            _bacterium.total_force = _bacterium.total_force + np.abs(
                                repulsion_x / len_pos * interaction_factor)

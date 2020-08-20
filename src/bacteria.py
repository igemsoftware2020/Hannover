#!/usr/bin/env python
# -*- coding: utf-8 -*-


# ********************************************************************************************
# imports
import random
import math
import numpy as np
import src.constants as c


# ********************************************************************************************
# main class <<<bacterium>>>
# the structure corresponds now to previous versions

class Bacterium(object):

    def __init__(self, position: np.ndarray = None,
                 width: float = c.BSUB_WIDTH,
                 length: float = c.BSUB_LENGTH,
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
            angle = [0, 0]
        else:
            self.angle = angle
        if position is None:
            self.position = [c.WINDOW_SIZE[0] / 2, c.WINDOW_SIZE / 2, 0]
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
        moment_of_inertia = c.BSUB_MASS / 6 * (3 * self.width ** 2 + self.length) + c.BSUB_MASS / 2 * self.width ** 2
        self.rotational_energy = 1/2 * moment_of_inertia * np.dot(self.velocity_angular, self.velocity_angular)
        return self.rotational_energy

    def get_translation_energy(self) -> float:
        return 1/2 * c.BSUB_MASS * np.dot(self.velocity, self.velocity)

    def grow(self):
        """
        grow bacteria for 1 second with speed growth_rate
        """
        # Make the bacteria grow
        # using a constant growth rate
        if self.living:
            self.length = self.length * (1 + c.BSUB_GROWTH_FACTOR)
            self.width = self.width * (1 + c.BSUB_GROWTH_FACTOR / 2)

    def random_cell_death(self):
        # Programmed cell death
        if random.random() > 1.0 - c.BSUB_MORTALITY_RATE:
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
            return angle + (0.5 - random.random()) * np.pi * 0.5

        def get_daughter_position(position, split_distance, angle):
            offset = (split_distance * math.sin(angle[0]) * math.cos(angle[1]),
                      split_distance * math.cos(angle[0]) * math.cos(angle[1]),
                      split_distance * math.sin(angle[1]))
            position = position + offset * volume_ratio
            return position

        volume_ratio = 0.4 + 0.2 * random.random()
        # Create daughter bacterium from self

        # Update parameters of daughter and mother bacterium
        daughter_bac_length = volume_ratio * self.length
        daughter_bac_angle = update_angle(self.angle)
        daughter_bac_position = get_daughter_position(position=self.position, split_distance=self.length*0.2,
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

    def get_volume(self):
        """ gives out the cubic volume equivalent """
        return self.width * self.width * self.length

    def move(self):
        """
        Move Bacteria for 1 time unit and add random movement and rotation
        :return:
        """
        # Set velocities depending on bacteria state
        # Active motion : ballistic movement
        if self.moving and self.living:
            self.velocity = self.velocity * np.asarray([math.sin(self.angle[0]), math.cos(self.angle[0]), 0],
                                                       dtype=np.float64) * c.TIME_STEP
            self.velocity_angular[0] = self.velocity_angular[0] + (0.5 - random.random()) * 0.1 * c.TIME_STEP
            self.velocity_angular[1] = self.velocity_angular[1] + (0.5 - random.random()) * 0.001 * c.TIME_STEP
        # Passive motion : random movement
        if not self.moving and self.living:
            self.velocity_angular[0] = self.velocity_angular[0] + (0.5 - random.random()) * 0.01 * c.TIME_STEP
            self.velocity_angular[1] = self.velocity_angular[1] + (0.5 - random.random()) * 0.001 * c.TIME_STEP

        # slight z-brownian random drift
        if not self.at_boundary():
            self.velocity[2] = self.velocity[2] + (0.5 - random.random()) * 0.1 * c.TIME_STEP
            # And gravity
            self.velocity[2] = self.velocity[2] - 0.981 * 0.5 * c.TIME_STEP
            self.position = self.position + self.velocity * c.TIME_STEP
            self.angle = self.angle + np.sqrt(np.dot(self.velocity_angular, self.velocity_angular)) * c.TIME_STEP
        elif self.at_boundary() == "X":
            self.velocity[0] = - self.velocity[0]
        elif self.at_boundary() == "Y":
            self.velocity[1] = - self.velocity[1]

    def at_boundary(self):
        x, y, z = self.position
        if x + self.length >= c.WINDOW_SIZE[0] or x - self.length <= 0:
            # elastic scattering at boundary
            return "X"
        elif y + self.length >= c.WINDOW_SIZE[1] or y - self.length <= 0:
            return "Y"
        return False

    def is_split_ready(self):
        def gaussian_distribution(x, mu=c.BSUB_CRITICAL_LENGTH, sigma=c.BSUB_CRITICAL_LENGTH * 0.12):
            return 1 / np.sqrt(2 * np.pi * sigma ** 2) * np.exp(- (x - mu) ** 2 / (2 * sigma ** 2))
        splitting_lengths = np.random.normal(c.BSUB_CRITICAL_LENGTH, c.BSUB_CRITICAL_LENGTH * 0.12)
        if splitting_lengths <= self.length <= splitting_lengths:
            probability_to_split = gaussian_distribution(self.length)
            return np.random.choice([True, False], p=[probability_to_split, 1-probability_to_split])
        else:
            return False

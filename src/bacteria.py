#!/usr/bin/env python
# -*- coding: utf-8 -*-


import math
# ********************************************************************************************
# imports
import random
from typing import Dict

import numpy as np

from src.constants import Constants as C
from src.utils import stokes_drag_force, gravitational_force, apply_rotation, rotation_matrix_y, rotation_matrix_x


# ********************************************************************************************
# main class <<<bacterium>>>

class Bacterium:

    def __init__(self, position: np.ndarray = None,
                 width: float = C.BSUB_WIDTH,
                 length: float = C.BSUB_LENGTH,
                 velocity: np.ndarray = np.asarray([0, 0, 0]),
                 angle: np.ndarray = None, force: np.ndarray = None,
                 living: bool = True, moving: bool = False, attached_to_surface: bool = False):
        """
        initialize a instance of the Bacteria class
        :param position: position of bacteria center [x_pox, y_pos, z_pos]
        :param width: width of ellipse in meter
        :param length:  length of ellipse in meter, default value 2 µm for B. subtilis
        :param velocity: velocity of bacteria [v_x, v_y, v_z] in m/s
        :param angle: angle of bacteria  measured to x axis in radian

        """

        if angle is None:
            self.angle = [random.randint(0, 360), random.randint(0, 360)]
        else:
            self.angle = angle

        self.position = position
        self.velocity = np.asarray([velocity[0], velocity[1], velocity[2]], dtype=np.int64)
        # rotate velocity in direction of orientation
        self.velocity = apply_rotation(self.velocity, rotation_matrix_x(self.angle[0]))
        self.velocity = apply_rotation(self.velocity, rotation_matrix_y(self.angle[1]))

        self.width = width
        self.length = length
        self.living = living
        self.moving = moving
        self.attached_to_surface = attached_to_surface

        self.velocity_angular = [0, 0, 0]

        self.force = force
        self.total_force = np.linalg.norm(self.force)

        self.rotational_energy = 0
        self.translation_energy = 0
        self.update_rotational_energy()
        self.update_translation_energy()
        self.total_energy = self.translation_energy + self.rotational_energy

    def update_rotational_energy(self):
        moment_of_inertia = C.BSUB_MASS / 6 * (3 * self.width ** 2 + self.length) + C.BSUB_MASS / 2 * self.width ** 2
        self.rotational_energy = 1 / 2 * moment_of_inertia * np.dot(self.velocity_angular, self.velocity_angular)

    def update_translation_energy(self):
        self.translation_energy = 1 / 2 * C.BSUB_MASS * np.dot(self.velocity, self.velocity)

    def get_volume(self):
        """ gives out the cubic volume equivalent """
        return self.width * self.width * self.length

    def grow(self):
        """
        grow bacteria for 1 second with speed growth_rate
        """
        # Make the bacteria grow
        # using a constant growth rate
        # TODO Make volume per time
        if self.living is True:
            self.length = self.length * (C.BSUB_GROWTH_FACTOR + 1)
        else:
            # if cell is dead, constant length
            pass

    def random_cell_death(self):
        # Programmed cell death
        if random.random() > 1.0 - C.BSUB_MORTALITY_RATE:
            self.living = False

    def __eq__(self, other):

        if get_bacteria_dict(self) == get_bacteria_dict(other):
            return True
        return False

    def split(self):
        """
        split bacteria  create new daughter bacteria with new values and update the mother bacterium
        :return: daughter bacteria
        """

        # Calculate new position and angle
        # Advanced new position: random angular component,
        #                       radial component sum of two radii*0.8

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
        daughter_bac_angle = self.angle  # same orientation?
        daughter_bac_position = get_daughter_position(position=self.position, split_distance=self.length * 0.2,
                                                      angle=daughter_bac_angle)
        daughter_bac = Bacterium(daughter_bac_position, self.width, daughter_bac_length,
                                 -self.velocity, self.angle, moving=True, force=-self.force)

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

    def update_velocity(self, dt=C.TIME_STEP):
        """
        Update velocity direction and value based on the acting force.
        Add Brownian movement in x,y,z direction
        Add random angle movement
        """
        if self.at_boundary() == 'X':
            apply_rotation(self.velocity, rotation_matrix_x(theta=np.pi))
        elif self.at_boundary() == 'Y':
            apply_rotation(self.velocity, rotation_matrix_y(theta=np.pi))

        acceleration = self.force / C.BSUB_MASS * 1E-6  # [um / s ** 2]
        # update velocities
        self.velocity[0] = self.velocity[0] + acceleration[0] * dt
        self.velocity[1] = self.velocity[1] + acceleration[1] * dt
        self.velocity[2] = self.velocity[2] + acceleration[2] * dt

        if self.position[2] < 1E-1:
            self.velocity[2] = 0
            self.position[2] = 0

        # self.velocity = apply_rotation(self.velocity, rotation_matrix_x(self.angle[0]))
        # self.velocity = apply_rotation(self.velocity, rotation_matrix_y(self.angle[1]))

        # add brownian movement
        self.velocity[0] = self.velocity[0] + random.uniform(-1, 1)  # add random velocity up to 1 um / s
        self.velocity[1] = self.velocity[1] + random.uniform(-1, 1)
        self.velocity[2] = self.velocity[2] + random.uniform(-1, 1)

        # update angular velocity
        # 3D  instantaneous angular velocity vector w = r x v / |r|^2
        self.velocity_angular = np.cross(self.position, self.velocity) / np.linalg.norm(self.position) ** 2
        # add random rotational velocity
        self.velocity_angular[0] = self.velocity_angular[0] + random.uniform(- 0.02, 0.02)
        self.velocity_angular[1] = self.velocity_angular[1] + random.uniform(- 0.02, 0.02)
        self.velocity_angular[2] = self.velocity_angular[2] + random.uniform(- 0.02, 0.02)

    def update_position(self, dt=C.TIME_STEP):
        """ update bacterium position based on velocity """
        self.position[0] = self.position[0] + self.velocity[0] * dt
        self.position[1] = self.position[1] + self.velocity[1] * dt
        if (self.position[2] + self.velocity[2] * dt) < 0:
            self.position[2] = 0
        else:
            self.position[2] = self.position[2] + self.velocity[2] * dt

        # update orientation
        self.angle[0] = self.angle[0] + self.velocity_angular[0] * dt
        self.angle[1] = self.angle[1] + self.velocity_angular[1] * dt
        self.angle[2] = self.angle[2] + self.velocity_angular[2] * dt

    def at_boundary(self):
        x, y, z = self.position
        if x + self.length >= C.WINDOW_SIZE[0] or x - self.length <= 0:
            # elastic scattering at boundary
            return "X"
        elif y + self.length >= C.WINDOW_SIZE[1] or y - self.length <= 0:
            return "Y"
        return False

    def is_split_ready(self):
        """
        checks if size of bacterium is long enough for splitting
        If bacterium is big enough, splitting occurs with a probability of 0.8
        returns True if splitting is possible, False otherwise
        """
        splitting_lengths = random.randrange(C.BSUB_CRITICAL_LENGTH - 1, C.BSUB_CRITICAL_LENGTH + 1)
        if splitting_lengths <= self.length:
            return np.random.choice([True, False], p=[0.8, 0.2])
        else:
            return False

    def update_acting_force(self):
        # Stokes drag force
        self.force = stokes_drag_force(radius=self.length, velocity=self.velocity)
        if (self.position[2] < 4) and self.attached_to_surface:
            # if distance from surface greater than 4 µm, add adhesion force
            self.force = self.force + C.MAX_CELL_SUBSTRATE_ADHESION * 1E-6 * np.asarray([0, 0, -1])
        # add gravitational force
        self.force += gravitational_force(C.BSUB_MASS)


def get_bacteria_dict(bacterium: Bacterium) -> Dict:
    """ returns the dict entry of a bacteria """
    return dict(
        position=[bacterium.position.tolist()],
        velocity=[bacterium.velocity.tolist()],
        angle=[bacterium.angle.tolist()],
        force=[bacterium.force.tolist()],
        total_force=[bacterium.total_force],
        total_energy=[bacterium.total_energy],
        living=[bacterium.living],
        moving=[bacterium.moving],
        length=[bacterium.length],
        width=[bacterium.width])

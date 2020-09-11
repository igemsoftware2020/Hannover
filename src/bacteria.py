#!/usr/bin/env python
# -*- coding: utf-8 -*-


import math
# ********************************************************************************************
# imports
import random
from typing import Dict

import numpy as np

from src.constants import Constants as C
from src.utils import stokes_drag_force


# ********************************************************************************************
# main class <<<bacterium>>>


class Bacterium:

    def __init__(self, position: np.ndarray = None,
                 width: float = C.BSUB_WIDTH,
                 length: float = C.BSUB_LENGTH,
                 velocity: np.ndarray = None,
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
            angle = [random.randint(0, 360), random.randint(0, 360)]
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
        self.attached_to_surface = attached_to_surface

        self.velocity_angular = [0, 0]
        self.force = force
        self.total_force = np.linalg.norm(self.force)

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
        # TODO Make volume per time
        if self.living is True:
            self.length = self.length * (C.BSUB_GROWTH_FACTOR + 1)
        else:
            # self.width  = self.width *gr_d_factor
            self.length = self.length * (1 - C.gr_d_factor)

    def random_cell_death(self):
        # Programmed cell death
        if random.random() > 1.0 - C.BSUB_MORTALITY_RATE:
            self.living = False

    def __eq__(self, other):

        if bacteria_dict(self) == bacteria_dict(other):
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
                                 self.velocity, self.angle, moving=True)

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
        if (self.at_boundary() == 'X') or (self.at_boundary() == 'Y'):
            # turn velocity by 90 degree
            return

        # projection of forces on movement direction
        # update velocities

        # add brown
        self.position[0] = self.position[0] + self.velocity[0] * dt
        self.position[1] = self.position[1] + self.velocity[1] * dt
        self.position[2] = self.position[2] + self.velocity[2] * dt
        self.angle[0] = self.angle[0] + self.velocity_angular[0] * dt
        self.angle[1] = self.angle[1] + self.velocity_angular[1] * dt

    def move(self, frameDim=C.WINDOW_SIZE, dt=C.TIME_STEP, friction=0.1):
        """
        Move Bacteria for 1 time unit and add random movement and rotation
        :return:
        """
        # Active motion
        if (self.moving == True) and (self.living == True):
            self.velocity[0] = self.velocity[0] + 1.0 * math.sin(self.angle[0]) * dt
            self.velocity[1] = self.velocity[1] + 1.0 * math.cos(self.angle[0]) * dt
            self.velocity_angular[0] = self.velocity_angular[0] + (0.5 - random.random()) * 0.1 * dt
            self.velocity_angular[1] = self.velocity_angular[1] + (0.5 - random.random()) * 0.001 * dt
        # Passive motion
        if (self.moving == False) and (self.living == True):
            self.velocity_angular[0] = self.velocity_angular[0] + (0.5 - random.random()) * 0.01 * dt
            self.velocity_angular[1] = self.velocity_angular[1] + (0.5 - random.random()) * 0.001 * dt

        # slight z-brownian random drift
        self.velocity[2] = self.velocity[2] + (0.5 - random.random()) * 0.1 * dt
        # And gravity
        self.velocity[2] = self.velocity[2] - 0.981 * 0.5 * dt

        # Boundary collision
        for position in self.get_position():
            if (False):
                xborder = frameDim[1] * 0.5
                if (position[0] < -xborder + self.width):
                    self.velocity[0] = self.velocity[0] + (-position[0] - xborder + self.width) * 0.1
                    self.total_force = self.total_force + np.linalg.norm(
                        -position[0] - xborder + self.width) * 0.1
                if (position[0] > xborder - self.width):
                    self.velocity[0] = self.velocity[0] + (-position[0] + xborder - self.width) * 0.1
                    self.total_force = self.total_force + np.linalg.norm(
                        -position[0] + xborder - self.width) * 0.1
            yborder = frameDim[0] * 0.5
            if (position[1] < -yborder + self.width):
                self.velocity[1] = self.velocity[1] + (-position[1] - yborder + self.width) * 0.1
                self.total_force = self.total_force + np.linalg.norm(
                    -position[1] - yborder + self.width) * 0.1
            #    if(position[1]>yborder-self.width):
            #        self.velocity[1] = self.velocity[1] + (-position[1]+yborder-self.width)*0.1
            #        self.totalForce_equivalent = self.totalForce_equivalent + absolute(-position[1]+yborder-self.width)*0.1
            # Bottom-boundary (z<=0)
            if (position[2] < -0 + self.width):
                self.velocity[2] = self.velocity[2] + (-position[2] + self.width) * 0.1
                # self.totalForce_equivalent = self.totalForce_equivalent + absolute(-position[2]+self.width)*0.1

                # Bounding-box torque
                positions = self.get_position()
                lenPos = len(positions)
                for index in range(lenPos - 1):
                    position = positions[index]
                    t_radius = (index - lenPos * 0.5)
                    # dx = _Bacterium.pos[0] - position[0]#self.pos[0]
                    # dy = _Bacterium.pos[1] - position[1]#self.pos[1]
                    dz = -position[2] + self.width  # self.pos[2]
                    repulsion_z = -dz  #
                    # dr = dx*dx+dy*dy+dz*dz
                    interactionfactor = 0.005
                    self.velocity_angular[1] = self.velocity_angular[1] - t_radius * math.cos(
                        self.angle[1]) * repulsion_z / lenPos * 0.05 * interactionfactor * dt

        self.position[0] = self.position[0] + self.velocity[0] * dt
        self.position[1] = self.position[1] + self.velocity[1] * dt
        self.position[2] = self.position[2] + self.velocity[2] * dt
        self.angle[0] = self.angle[0] + self.velocity_angular[0] * dt
        self.angle[1] = self.angle[1] + self.velocity_angular[1] * dt

        self.velocity[0] = self.velocity[0] * friction
        self.velocity[1] = self.velocity[1] * friction
        self.velocity[2] = self.velocity[2] * friction
        self.velocity_angular[0] = self.velocity_angular[0] * friction
        self.velocity_angular[1] = self.velocity_angular[1] * friction

        # self.velocity_angular[1] = self.velocity_angular[1]+0.1
        # Total Force sum
        # kind of proportional to biofilm pressure
        # decays over time
        self.total_force = self.total_energy * 0.9

    def at_boundary(self):
        x, y, z = self.position
        if x + self.length >= C.WINDOW_SIZE[0] or x - self.length <= 0:
            # elastic scattering at boundary
            return "X"
        elif y + self.length >= C.WINDOW_SIZE[1] or y - self.length <= 0:
            return "Y"
        return False

    def is_split_ready(self):

        splitting_lengths = random.randrange(C.BSUB_CRITICAL_LENGTH - 1, C.BSUB_CRITICAL_LENGTH + 1)
        if splitting_lengths <= self.length:
            return np.random.choice([True, False], p=[0.8, 0.2])
        else:
            return False

    def update_acting_force(self):
        # Stokes drag force
        self.force = stokes_drag_force(radius=self.length, velocity=self.velocity)
        if (self.position[3] < 4) and self.attached_to_surface:
            # if distance from surface greater than 4 µm, add adhesion force
            self.force = self.force + C.MAX_CELL_SUBSTRATE_ADHESION


def bacteria_dict(bacterium: Bacterium) -> Dict:
    """ returns the dict entry of a bacteria """
    return dict(
        position=[bacterium.position.tolist()],
        velocity=[bacterium.velocity.tolist()],
        angle=[bacterium.angle],
        force=[bacterium.force],
        total_force=[bacterium.total_force],
        total_energy=[bacterium.total_energy],
        living=[bacterium.living],
        moving=[bacterium.moving],
        length=[bacterium.length],
        width=[bacterium.width])

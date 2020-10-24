#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ********************************************************************************************
# imports
import math
import random
from typing import Dict

# custom libraries
import BiofilmSimulation.constants as c
import numpy as np
import scipy.stats
from BiofilmSimulation.formulas import stokes_drag_force, gravitational_force, apply_rotation, rotation_matrix_y, \
    rotation_matrix_x, \
    lennard_jones_force, rotation_matrix_z


# ********************************************************************************************
# main class <<<bacterium>>>

class Bacterium:

    def __init__(self, constants: c, strain: str = None, position: np.ndarray = None,
                 velocity: np.ndarray = np.asarray([np.random.normal(0, 0.5),
                                                    np.random.normal(0, 0.5),
                                                    np.random.normal(0, 0.5)]),
                 angle: np.ndarray = None, force: np.ndarray = None,
                 living: bool = True, moving: bool = False,
                 attached_to_surface: bool = False, length: float = None):
        """
        initialize a instance of the Bacteria class
        :param constants: c used for the bacterium. Object of c class
        :param strain:  str can be set to "E.Coli" or "B.Sub.". Default uses Bacteria type selected in constants
        :param position: position of bacteria center [x_pox, y_pos, z_pos]
        :param velocity: velocity of bacteria [v_x, v_y, v_z] in m/s
        :param angle: angle of bacteria  measured to x axis in degree
        :param force: acting force of bacteria in each direction in N
        :param living: True if bacteria is alive, false else
        :param moving: True if bacteria is moving, false else
        :param attached_to_surface: True if bacteria is attached to surface, false else
        :param length:  length of ellipse in meter, default value 2 µm for B. sub
        """
        self.constants = constants
        # initial position
        self.position = position
        # have to add this here, so it will be stored in the log file
        self.height = self.position[2]
        self.velocity: np.ndarray = np.asarray([velocity[0], velocity[1], velocity[2]], dtype=np.int64)

        # initial orientation
        if angle is None:
            self.angle = [random.randint(0, 360), random.randint(0, 360)]
        else:
            self.angle = angle

        # rotate velocity in direction of orientation
        self.velocity: np.ndarray = apply_rotation(self.velocity, rotation_matrix_x(self.angle[0]))
        self.velocity: np.ndarray = apply_rotation(self.velocity, rotation_matrix_y(self.angle[1]))

        self.mass = self.constants.get_bsub_constants(key="MASS")
        self.length = length

        if strain is None:
            self.strain = constants.bac_type
        else:
            self.strain = strain
        # set remaining parameters according to selected bacteria strain
        if self.strain == "B.Sub.":
            self.width = self.constants.get_bsub_constants(key="WIDTH")
            if length is None:
                self.length = self.constants.get_bsub_constants(key="LENGTH")

            self.growth_rate = self.constants.get_bsub_constants(key="GROWTH_RATE")
            self.mortality_rate = self.constants.get_bsub_constants(key="MORTALITY_RATE")
            self.critical_length = self.constants.get_bsub_constants(key="CRITICAL_LENGTH")
        elif self.strain == "E.Coli.":
            self.width = self.constants.get_ecoli_constants(key="WIDTH")
            if length is None:
                self.length = self.constants.get_ecoli_constants(key="LENGTH")

            self.growth_rate = self.constants.get_ecoli_constants(key="GROWTH_RATE")
            self.mortality_rate = self.constants.get_ecoli_constants(key="MORTALITY_RATE")
            self.critical_length = self.constants.get_ecoli_constants(key="CRITICAL_LENGTH")

        self.update_mass()
        self.living = living
        self.moving = moving
        self.attached_to_surface = attached_to_surface

        self.velocity_angular = [0, 0, 0]

        self.force: np.ndarray = force
        self.total_force = np.linalg.norm(self.force)
        self.acceleration = self.force / self.mass * 1E-6

        self.rotational_energy = 0
        self.translation_energy = 0
        self.update_rotational_energy()
        self.update_translation_energy()
        self.total_energy = self.translation_energy + self.rotational_energy

    def __eq__(self, other):
        """ new __eq__ based on bacteria parameters """
        if get_bacteria_dict(self) == get_bacteria_dict(other):
            return True
        return False

    def update_velocity(self):
        """
        Update velocity direction and value based on the acting force.
        Add Brownian movement in x,y,z direction
        Add random angle movement
        """
        dt = self.constants.get_simulation_constants(key="time_step")
        # update velocities
        self.velocity[0] += self.acceleration[0] * dt
        self.velocity[1] += self.acceleration[1] * dt
        self.velocity[2] += self.acceleration[2] * dt

        # rotate velocity in direction of orientation
        # self.velocity: np.ndarray = apply_rotation(self.velocity, rotation_matrix_x(self.angle[0]))
        # self.velocity: np.ndarray = apply_rotation(self.velocity, rotation_matrix_y(self.angle[1]))
        # self.velocity: np.ndarray = apply_rotation(self.velocity, rotation_matrix_z(self.angle[2]))

    def update_position(self):
        """ update bacterium position based on velocity """
        dt = self.constants.get_simulation_constants(key="time_step")
        self.position[0] += self.velocity[0] * dt + 1 / 2 * self.acceleration[0] * dt ** 2
        self.position[1] += self.velocity[1] * dt + 1 / 2 * self.acceleration[1] * dt ** 2
        self.position[2] += self.velocity[2] * dt + 1 / 2 * self.acceleration[2] * dt ** 2

        local_rnd_1 = np.random.RandomState()
        local_rnd_2 = np.random.RandomState()
        local_rnd_3 = np.random.RandomState()
        self.position[0] = local_rnd_1.normal(loc=self.position[0], scale=0.5)
        self.position[1] = local_rnd_2.normal(loc=self.position[1], scale=0.5)
        self.position[2] = local_rnd_3.normal(loc=self.position[2], scale=0.01)

        if self.position[2] < self.length:
            self.position[2] = self.width
            self.attached_to_surface = True

    def update_orientation(self):
        """ update bacterium orientation """
        # update angular velocity
        # 3D  instantaneous angular velocity vector w = r x v / |r|^2
        self.velocity_angular = np.cross(self.position, self.velocity) / np.linalg.norm(self.position) ** 2
        # add random rotational velocity
        local_rnd_1 = np.random.RandomState()
        local_rnd_2 = np.random.RandomState()
        local_rnd_3 = np.random.RandomState()

        self.velocity_angular[0] += local_rnd_1.normal(loc=self.velocity_angular[0], scale=0.5)
        self.velocity_angular[1] += local_rnd_2.normal(loc=self.velocity_angular[1], scale=0.5)
        self.velocity_angular[2] += local_rnd_3.normal(loc=self.velocity_angular[2], scale=0.5)

        local_rnd_1 = np.random.RandomState()
        local_rnd_2 = np.random.RandomState()
        local_rnd_3 = np.random.RandomState()

        self.angle[0] = local_rnd_1.normal(loc=self.angle[0], scale=0.2) + self.velocity_angular[0]
        self.angle[1] = local_rnd_2.normal(loc=self.angle[1], scale=0.2) + self.velocity_angular[1]
        self.angle[2] = local_rnd_3.normal(loc=self.angle[2], scale=0.2) + self.velocity_angular[2]

    def update_acting_force(self):
        """
        Calculates all forces acting on the bacteria
        and updates the according parameter of the bacteria.
        Forces included:
         Stokes drag force, bacterium- bacterium adhesion,
        bacterium-Substrate adhesion, gravitation
        """
        # offset
        # self.force = self.mass * self.acceleration * 1E6
        # Stokes drag force
        self.force = 0
        self.force = np.add(self.force, stokes_drag_force(radius=self.length, velocity=self.velocity,
                                                          viscosity=self.constants.EFFECTIVE_VISCOSITY_EPS)
                            )
        # add gravitational force
        self.force = np.add(self.force, gravitational_force(self.mass))
        self.force = np.add(self.force, bac_substrate_interaction_force(self))
        self.total_force = np.linalg.norm(self.force)

    def update_acceleration(self):
        """ calculates and sets acceleration """
        self.acceleration = self.force / self.mass * 1E-6

    def update_rotational_energy(self):
        """ updates the rotational energy """
        moment_of_inertia = self.mass / 6 * (3 * self.width ** 2 + self.length) + self.mass / 2 * self.width ** 2
        self.rotational_energy = 1 / 2 * moment_of_inertia * np.dot(self.velocity_angular, self.velocity_angular)

    def update_mass(self):
        """update mass of bacteria on based on volume"""
        # Mass / Volume ration for grown bacteria
        ratio = self.constants.get_bac_constants(key="MASS") / (
                np.pi * self.constants.get_bac_constants(key="WIDTH") ** 2
                * self.constants.get_bac_constants(key="LENGTH"))
        volume = np.pi * self.length * self.width ** 2
        self.mass = ratio * volume

    def update_translation_energy(self):
        """ updates the translation energy """
        self.translation_energy = 1 / 2 * self.mass * np.dot(self.velocity, self.velocity)

    def split(self):
        """
        split bacteria  create new daughter bacteria with new values and update the mother bacterium
        :return: daughter bacteria
        """

        # Calculate new position and angle
        # Advanced new position: random angular component,
        #                       radial component sum of two radii*0.8

        def get_daughter_position(mother_bac: Bacterium):
            r_mother = mother_bac.position
            v_mother = mother_bac.velocity
            split_distance = mother_bac.length / 2 + 0.5  # lengths in microns
            r_daughter = r_mother + v_mother / np.linalg.norm(v_mother) * split_distance
            return r_daughter

        # Create daughter bacterium from self
        # Update parameters of daughter and mother bacterium
        # make splitting ration from normal distrbiturio
        daughter_bac_position = get_daughter_position(self)
        daughter_bac_position[2] = self.position[2]
        daughter_bac_length = self.length * np.random.normal(0.5, 0.5 * 0.07)
        daughter_bac_velocity = np.asarray(self.velocity / 2)
        daughter_bac_force = np.asarray(self.force / 2)

        daughter_bac = Bacterium(constants=self.constants, strain=self.strain,
                                 angle=self.angle, force=daughter_bac_force,
                                 living=True, moving=True,
                                 attached_to_surface=self.attached_to_surface,
                                 velocity=daughter_bac_velocity,
                                 position=daughter_bac_position, length=daughter_bac_length)

        daughter_bac.acceleration = self.acceleration / 2
        daughter_bac.update_mass()
        daughter_bac.update_acting_force()
        daughter_bac.update_acceleration()
        daughter_bac.update_velocity()
        daughter_bac.update_position()

        # update mother cell
        self.length = self.length  - daughter_bac_length
        self.velocity = self.velocity / 2
        self.force = self.force / 2
        self.update_mass()
        return daughter_bac

    def get_volume(self):
        """ gives out the cubic volume equivalent """
        return self.width * self.width * self.length

    def grow(self):
        """
        grow bacteria for 1 second with speed growth_rate
        """
        # Make the bacteria grow
        # using a constant growth rate
        if self.living is True:
            self.length = self.length + self.growth_rate * self.constants.get_simulation_constants(key="time_step")
        else:
            # if cell is dead, constant length
            pass

    def random_cell_death(self):
        """ random cell dying """
        if random.random() > 1.0 - self.mortality_rate:
            self.living = False
            self.moving = False

    def detach(self):
        if (self.attached_to_surface is True) & (np.random.random() > 0.9):
            self.attached_to_surface = False
            self.acceleration[2] = 0.01

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

    def at_boundary(self):
        """ checks if bacteria is at the edge of the simulation plane"""
        x, y, z = self.position
        window_size = self.constants.get_simulation_constants(key="window_size")
        if x + self.length >= window_size[0] or x - self.length <= 0:
            # elastic scattering at boundary
            return "X"
        elif y + self.length >= window_size[1] or y - self.length <= 0:
            return "Y"
        return False

    def is_split_ready(self):
        """
        checks if size of bacterium is long enough for splitting
        If bacterium is big enough, splitting occurs with a probability
        of a normal distribution with mean at self.critical_length
        returns True if splitting is possible, False otherwise
        """
        probability = scipy.stats.norm.cdf(self.length, loc=self.critical_length, scale=self.critical_length * 0.12)
        return np.random.choice([True, False], p=[probability, 1 - probability])


def get_bacteria_dict(bacterium: Bacterium) -> Dict:
    """ returns the dict entry of a bacteria """
    return dict(
        position=[bacterium.position.tolist()],
        velocity=[bacterium.velocity.tolist()],
        acceleration=[bacterium.acceleration.tolist()],
        height=[bacterium.height],
        angle=[bacterium.angle.tolist()],
        force=[bacterium.force.tolist()],
        total_force=[bacterium.total_force],
        total_energy=[bacterium.total_energy],
        living=[bacterium.living],
        moving=[bacterium.moving],
        length=[bacterium.length],
        width=[bacterium.width],
        mass=[bacterium.mass]
    )


def bac_substrate_interaction_force(self: Bacterium):
    """
        returns force vector of bacterium substrate interaction
        """
    if self.position[2] > 4:
        force = lennard_jones_force(self.position[2], f_min=-self.constants.MAX_CELL_SUBSTRATE_ADHESION,
                                    r_min=1) \
                * np.asarray([0, 0, 1])
    else:
        force = self.constants.MAX_CELL_SUBSTRATE_ADHESION * np.asarray([0, 0, -1])
    return force


def bac_bac_interaction_force(self: Bacterium, other: Bacterium):
    """
        returns force vector of cell- cell interaction.
        Force direction is in direction of the distance vector between the bacteria.
        Force value based on Lennard-Jones Potential / Soft-repulsive potential
        """

    if np.linalg.norm(distance_vector(self, other)) > 1.9:
        distance_abs = np.linalg.norm(distance_vector(self, other))
        return distance_vector(self, other) / distance_abs * \
               lennard_jones_force(distance_abs, f_min=-self.constants.MAX_CELL_CELL_ADHESION, r_min=2)
    return self.constants.MAX_CELL_CELL_ADHESION * distance_vector(self, other) \
           / np.linalg.norm(distance_vector(self, other))


def distance_vector(self: Bacterium, other: Bacterium):
    """ return distance vector between two bacteria """
    return self.position - other.position

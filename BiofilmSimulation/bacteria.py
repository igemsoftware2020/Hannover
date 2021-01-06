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

    def __init__(self, constants: c,
                 strain: str = None,
                 position: np.ndarray = None,
                 velocity: np.ndarray = np.empty((3, 1)),
                 angle: np.ndarray = np.empty((3, 1)),
                 force: np.ndarray = np.empty((3, 1)),
                 living: bool = True,
                 attached_to_surface: bool = False,
                 length: float = None,
                 index: int = None):
        """
        initialize a instance of the Bacteria class
        :param constants: c used for the bacterium. Object of c class
        :param strain:  str can be set to "E.Coli" or "B.Sub.". Default uses Bacteria type selected in constants
        :param position: position of bacteria center [x_pox, y_pos, z_pos]
        :param velocity: velocity of bacteria [v_x, v_y, v_z] in m/s
        :param angle: angle of bacteria. Measured in xy plane and in zy-plane. Measured from x- axis in degree.
        :param force: acting force of bacteria in each direction in N
        :param living: True if bacteria is alive, false else
        :param attached_to_surface: True if bacteria is attached to surface, false else
        :param length:  length of ellipse in micrometer, default value 2 µm for B. sub
        """
        self.index: int = index
        self.constants: c = constants
        # initial position
        self.position: np.ndarray = position
        # have to add this here, so it will be stored in the log file
        self.velocity: np.ndarray = velocity

        # initial orientation
        self.angle = angle
        self.length = length

        # set remaining parameters according to selected bacteria strain
        self.width = self.constants.get_bac_constants(key="WIDTH")
        self.mass = self.constants.get_bac_constants(key="MASS")
        self.growth_rate = self.constants.get_bac_constants(key="GROWTH_RATE")
        self.mortality_rate = self.constants.get_bac_constants(key="MORTALITY_RATE")
        self.critical_length = self.constants.get_bac_constants(key="CRITICAL_LENGTH")

        if length is None:
            self.length = self.constants.get_bac_constants(key="LENGTH")

        if strain is None:
            self.strain = constants.bac_type
        else:
            self.strain = strain

        self.update_mass()
        self.living = living
        self.attached_to_surface = attached_to_surface

        self.velocity_angular = [0, 0, 0]  # [omega_x, omega_y, omega_z]

        self.force: np.ndarray = force  # [N]
        self.total_force = np.linalg.norm(self.force)  # [N]
        self.acceleration = self.force / self.mass * 1E-6  # [um / s^2]

        # Use this class parameters to check for energy conservation after splitting process
        self.rotational_energy = 0  # E_rot = 1/2 * I * omega^2
        self.translation_energy = 0  # E_kin = 1/2 * m * v^2
        self.total_energy = self.translation_energy + self.rotational_energy

    def __eq__(self, other):
        """ new __eq__ based on bacteria parameters """
        if self.index == other.index:
            return True
        return False

    def get_index(self):
        return self.index

    def update_velocity(self):
        """
        Update velocity direction and value based on the acting force.
        """
        dt = self.constants.get_simulation_constants(key="time_step")
        # update velocities
        self.velocity += self.acceleration * dt

        # rotate velocity in direction of orientation
        self.velocity: np.ndarray = apply_rotation(self.velocity, rotation_matrix_x(self.angle[0]))
        self.velocity: np.ndarray = apply_rotation(self.velocity, rotation_matrix_y(self.angle[1]))

    def update_position(self):
        """
        update bacterium position based on velocity,
        add brownian movement to position
         """
        dt = self.constants.get_simulation_constants(key="time_step")
        self.position += self.velocity * dt + 1 / 2 * self.acceleration * dt ** 2

        local_rnd_1 = np.random.RandomState()
        local_rnd_2 = np.random.RandomState()
        local_rnd_3 = np.random.RandomState()

        self.position[0] = local_rnd_1.normal(loc=self.position[0], scale=0.5)
        self.position[1] = local_rnd_2.normal(loc=self.position[1], scale=0.5)
        self.position[2] = local_rnd_3.normal(loc=self.position[2], scale=0.5)

        if self.position[2] < self.length:
            self.attached_to_surface = True

    def update_orientation(self):
        """
        update bacterium orientation,
        add brownian movement
         """
        # update angular velocity
        # 3D  instantaneous angular velocity vector w = r x v / |r|^2
        # TODO: Transformation in spherical and cartesian coordinates and back!
        self.velocity_angular = np.cross(self.position, self.velocity) / np.linalg.norm(self.position) ** 2

        dt = self.constants.get_simulation_constants(key="time_step")
        self.angle[0] += self.velocity_angular[0] * dt

        local_rnd_1 = np.random.RandomState()
        local_rnd_2 = np.random.RandomState()

        self.angle[0] = local_rnd_1.normal(loc=self.angle[0], scale=0.2) + self.velocity_angular[0]
        self.angle[1] = local_rnd_2.normal(loc=self.angle[1], scale=0.2) + self.velocity_angular[1]

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
        self.force = np.asarray([0, 0, 0])
        self.force = np.add(self.force, stokes_drag_force(radius=self.width / 2, velocity=self.velocity,
                                                          viscosity=self.constants.EFFECTIVE_VISCOSITY_EPS)
                            )
        # add gravitational force
        self.force = np.add(self.force, gravitational_force(self.mass))
        self.force = np.add(self.force, bac_substrate_interaction_force(self))
        self.total_force = np.linalg.norm(self.force)

    def update_acceleration(self):
        """ calculates and sets acceleration in [um / s^2]"""
        self.acceleration = self.force / self.mass * 1E6

    def update_rotational_energy(self):
        """ updates the rotational energy """
        moment_of_inertia = self.mass / 6 * 3 * (self.width / 2) ** 2 + (self.length / 2) + self.mass / 2 * (
                    self.width / 2) ** 2
        self.rotational_energy = 1 / 2 * moment_of_inertia * np.dot(self.velocity_angular, self.velocity_angular)

    def update_mass(self):
        """update mass of bacteria on based on volume"""
        # Mass / Volume ration for grown bacteria
        ratio = self.constants.get_bac_constants(key="MASS") / (
                np.pi * self.constants.get_bac_constants(key="WIDTH") ** 2
                * self.constants.get_bac_constants(key="LENGTH"))
        volume = np.pi * self.length * (self.width / 2) ** 2
        self.mass = ratio * volume

    def update_translation_energy(self):
        """ updates the translation energy """
        self.translation_energy = 1 / 2 * self.mass * np.dot(self.velocity, self.velocity) * 1E-12

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
        daughter_bac_length = self.length * np.random.normal(0.5, 0.5 * 0.07)
        daughter_bac_velocity = np.asarray(self.velocity / 2)

        daughter_bac = Bacterium(constants=self.constants, strain=self.strain,
                                 angle=self.angle,
                                 living=True,
                                 attached_to_surface=self.attached_to_surface,
                                 velocity=daughter_bac_velocity,
                                 position=daughter_bac_position, length=daughter_bac_length)

        daughter_bac.acceleration = self.acceleration / 2
        daughter_bac.update_mass()
        daughter_bac.update_acting_force()
        daughter_bac.update_acceleration()

        # update mother cell
        self.length = self.length - daughter_bac_length
        self.velocity = self.velocity / 2
        self.update_mass()
        self.update_acting_force()
        self.update_acceleration()

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
            self.length += self.growth_rate * self.constants.get_simulation_constants(key="time_step")
        else:
            # if cell is dead, constant length
            pass

    def random_cell_death(self):
        """ random cell dying """
        if random.random() > 1.0 - self.mortality_rate:
            self.living = False

    def detach(self):
        """ detach bacteria from surface. Models bacteria diffusing from substrate into space."""
        if (self.attached_to_surface is True) & (np.random.random() > 0.9):
            self.attached_to_surface = False
            self.acceleration[2] = 0.01

    def get_position(self) -> np.ndarray:
        """ discretize the bacteria position in space  """
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


# Functions, which depend on the Bacteria class


def get_bacteria_dict(bacterium: Bacterium) -> Dict:
    """ returns the dict entry of a bacteria """
    return dict(
        position=[bacterium.position.tolist()],
        velocity=[bacterium.velocity.tolist()],
        acceleration=[bacterium.acceleration.tolist()],
        angle=[bacterium.angle.tolist()],
        force=[bacterium.force.tolist()],
        total_force=[bacterium.total_force],
        total_energy=[bacterium.total_energy],
        living=[bacterium.living],
        length=[bacterium.length],
        width=[bacterium.width],
        mass=[bacterium.mass]
    )


def bac_substrate_interaction_force(self: Bacterium):
    """
        returns force vector of bacterium substrate interaction
        """
    if not self.attached_to_surface:
        force = lennard_jones_force(self.position[2], f_min=-self.constants.MAX_CELL_SUBSTRATE_ADHESION,
                                    r_min=self.length / 2) \
                * np.asarray([0, 0, -1])
    else:
        force = self.constants.MAX_CELL_SUBSTRATE_ADHESION * np.asarray([0, 0, -1])
    return force


def bac_bac_interaction_force(self: Bacterium, other: Bacterium):
    """
        returns force vector of cell- cell interaction.
        Force direction is in direction of the distance vector between the bacteria.
        Force value based on Lennard-Jones Potential / Soft-repulsive potential
        """

    distance_absolute = np.linalg.norm(distance_vector(self, other))

    if distance_absolute > 1.9:
        return distance_vector(self, other) / distance_absolute * \
               lennard_jones_force(distance_absolute, f_min=-self.constants.MAX_CELL_CELL_ADHESION, r_min=2)
    return self.constants.MAX_CELL_CELL_ADHESION * distance_vector(self, other) / distance_absolute


def distance_vector(self: Bacterium, other: Bacterium):
    """ return distance vector between two bacteria """
    return self.position - other.position

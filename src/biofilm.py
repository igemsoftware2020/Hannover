#!/usr/bin/env python
# -*- coding: utf-8 -*-


# ********************************************************************************************
# imports
import numpy as np
import random

# custom libraries
from src.bacteria import Bacterium
import src.constants as c


class Biofilm(object):

    def __init__(self):
        self.bacteria = []

    def spawn(self):
        for _ in range(c.START_NUMBER_BACTERIA):
            pos = np.asarray(
                [(random.random() - 0.5) * c.WINDOW_SIZE[0], (random.random() - 0.5) * c.WINDOW_SIZE[1], 0])
            velocity = np.asarray([random.gauss(c.BSUB_MEAN_SPEED, 10), random.gauss(c.BSUB_MEAN_SPEED, 10), 0])
            bac = Bacterium(position=pos, velocity=velocity)
            self.bacteria.append(bac)

    @staticmethod
    def check_energy_conservation(bacterium1: Bacterium, bacterium2: Bacterium, total_energy_before):
        if bacterium1.total_energy + bacterium2.total_energy != total_energy_before:
            raise ValueError("Energy conversation broken while splitting.")

    def evolve(self):
        for bacterium in self.bacteria:
            bacterium.grow()

            if not bacterium.living and bacterium.getVolume() < 3300:
                self.bacteria.remove(bacterium)
            # Manage repulsion
            # (The most expensive task)
            for _bacterium in self.bacteria:
                [_bacterium, bacterium] = Biofilm.interaction(bacterium, _bacterium)

            bacterium.move()

            # Manage Bacterial splitting
            # Add a little bit of random until -----it looks good and real-----
            if bacterium.is_split_ready() and bacterium.living:
                energy_before = bacterium.total_energy
                daughter = bacterium.split()
                self.check_energy_conservation(bacterium, daughter, energy_before)
                self.bacteria.append(daughter)

            # Manage Bacterial Motion-Mode
            # Enable/Disable motion mode
            if bacterium.living:
                if random.random() > 1.0 - c.MOTION_ACTIVATION_PROBABILITY:
                    bacterium.moving = True
                if random.random() > 1.0 - c.MOTION_DEACTIVATION_PROBABILITY:
                    bacterium.moving = False
            else:
                bacterium.moving = False

    @staticmethod
    def interaction(bacterium1: Bacterium, bacterium2: Bacterium):
        """ return interaction term with bacteria in local environment"""

        def lennard_jones_force(r, epsilon, r_min):
            return - epsilon * (12 * (r_min ** 12 / r ** 13) - 12 * (r_min ** 6 / r ** 7))

        distance_vector = bacterium1.position - bacterium2.position
        distance = np.sqrt(np.dot(distance_vector, distance_vector))

        if distance < bacterium1.length * 4 and distance != 0:
            # Values for epsilon and r_min kinda random
            repulsive_force = lennard_jones_force(distance, epsilon=0.1, r_min=bacterium1.length)
            acceleration = repulsive_force / c.BSUB_MASS
            def new_position(a, v, t, r): return 1 / 2 * a * t ** 2 + v * t + r
            bacterium1.position = new_position(acceleration, bacterium1.velocity, c.TIME_STEP, bacterium1.position)
            bacterium2.position = new_position(acceleration, bacterium2.velocity, c.TIME_STEP, bacterium2.position)

        return bacterium1, bacterium2

    def sort_by_depth(self, axis, _reverse):
        sorted_bacteria = self.bacteria
        # To return a new list, use the sorted() built-in function...
        return sorted(sorted_bacteria, key=lambda x: x.position[axis], reverse=_reverse)

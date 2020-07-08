#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a awesome
        python script!
"""

import numpy as np
from matplotlib import pyplot as plt


class Bacteria(object):
    growth_rate = 1.3  # um / s , use Monod equation if differenz growth rates with substrate concentration
    split_length = 20   # cell splits, if height is greater that 10 um
    aspect_ratio = 0.4

    def __init__(self, pos, length, velocity, angle, index):
        """
        initialize a instance of the Bacteria class
        :param pos: position of bacteria center [x_pox, y_pos, z_pos]
        :param width: diameter of ellipse
        :param length:  length of ellipse
        :param velocity: velocity of bacteria [v_x, v_y, v_z]
        :param angle: angle of bacteria in radian  measured to x axis
        :param index: integer value, index of bacteria
        """

        self.pos = pos
        self.length = length  # lenghts in um
        self.width = length * self.aspect_ratio    # bacteria is a ellipse with radii width / length
        self.angle = angle    # angle is orientation of cell measured to x- axes in degree
        self.velocity = velocity[0], velocity[1]
        self.index = index

    def grow(self):
        """
        grow bacteria for 1 second with speed growth_rate
        """
        self.length = self.length * self.growth_rate
        self.width = self.length * self.aspect_ratio

    def split(self):
        """
        split bacteria if critical length is reached, create new bacteria with new values
        :return: daughter bacteria
        """

        def new_centers(center, angle, distance):
            """
            calculate new center of bacteria after splitting, calculation take into account,
            that splitting will more often appear along long axis
             """
            def rotate_vector(vector, angle):
                rotation_matrix = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
                rotated = np.dot(rotation_matrix, vector)
                return rotated / np.linalg.norm(rotated)

            angle_to_y = np.arccos(np.dot(center, np.array([0, 1])) / (np.linalg.norm(center)))
            if center[0] < 0:
                angle_to_y = -angle_to_y
            center1 = center
            center2 = center + distance * rotate_vector(rotate_vector(center, angle), angle_to_y)
            return center1, center2

        # update parameters of mother
        self.length = self.length / 2
        self.width = self.width / 2
        self.velocity = self.velocity[0] / 2, self.velocity[1] / 2

        # calculate new centers
        c1, c2 = new_centers(self.pos, self.angle, 2 * self.length)
        # set new bacteria
        bac_1 = Bacteria(pos=c1,  length=self.length, velocity=self.velocity,
                         angle=self.angle, index=self.index)
        bac_2 = Bacteria(pos=c2, length=self.length,
                         velocity=[-self.velocity[0], -self.velocity[1]],
                         angle=self.angle, index=self.index + 1)
        del self
        return bac_1, bac_2

    def move(self):
        """
        Move Bacteria for 1 time unit and add random movement and rotation
        :return:
        """
        rand_x, rand_y = self.brownian(100)
        self.pos = np.array([self.pos[0] + self.velocity[0] * 1 + rand_x[0],
                             self.pos[1] + self.velocity[1] * 1 + rand_y[0]])

        self.angle = self.angle * np.random.normal(loc=1, scale=2, size=1)[0]

    def interaction(self):
        """ return interaction term with bacteria in local environment"""

    def print_info(self):
        """
        print all values of parameters of bacteria
        :return:
        """
        print("\n**** Current values of bacci {index} ****\n"
              "  Position       ({x_pos:f}, {y_pos:f})\n"
              "  Width, Length      {width:.2f},{length:.2f}\n"
              "  angle            {o:.1f}\n"
              "  velocity           {vx:.2f}, {vy:.2f}".format(index=self.index, x_pos=self.pos[0],
                                                               y_pos=self.pos[1],
                                                               width=self.width, length=self.length,
                                                               o=self.angle, vx=self.velocity[0], vy=self.velocity[1])
              )

    @staticmethod
    def brownian(n):
        """ simulating random movement for n steps, returns 2 times n random numbers"""
        x = np.cumsum(np.random.randn(n))
        y = np.cumsum(np.random.randn(n))
        return x, y

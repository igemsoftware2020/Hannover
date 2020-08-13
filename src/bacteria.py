#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a awesome
        python script!
"""
# ********************************************************************************************
# imports
import random
import math


# ********************************************************************************************
# general dependency functions
def absolute(value):
    if (value >= 0):
        return value
    else:
        return -value


# ********************************************************************************************
# main class <<<bacterium>>>
# the structure corresponds now to previous versions

class bacterium(object):

    def __init__(self, pos, width, length, velocity, angle, totalForce_equivalent, living, moving):
        """
        initialize a instance of the Bacteria class
        :param pos: position of bacteria center [x_pox, y_pos, z_pos]
        :param width: width of ellipse
        :param length:  length of ellipse
        :param velocity: velocity of bacteria [v_x, v_y, v_z]
        :param angle: angle of bacteria in radian  measured to x axis
        """

        self.pos = [pos[0], pos[1], pos[2]]

        self.width = width  # bacteria is a ellipse with radii width , length, width
        self.length = length  # lenghts in um

        self.angle = angle  # angle is orientation of cell measured to x- axes in degree
        self.velocity_angular = [0, 0]
        self.velocity = [velocity[0], velocity[1], velocity[2]]

        self.totalForce_equivalent = totalForce_equivalent
        self.living = living
        self.moving = moving

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
                                                               o=self.angle, vx=self.v_x, vy=self.v_y)
              )

    def grow(self, gr_factor, gr_pr_i, gr_factor_inv, gr_d_factor, mortility_rate):
        """
        grow bacteria for 1 second with speed growth_rate
        """
        # Make the bacteria grow
        # using a constant growth rate ( volume per time)
        # growths inhibition factor
        if (self.living == True):
            if (self.moving == False):
                growth_suppressor = gr_factor * gr_pr_i / (gr_pr_i + self.totalForce_equivalent * 0.5) - gr_factor_inv
                volume = self.getVolume()
                growth_factor = (volume + volume ** (1 / 3) * 5.0 * growth_suppressor * growth_suppressor) / volume
                # self.width  = self.width *growth_factor
                self.length = self.length * growth_factor
        else:
            # self.width  = self.width *gr_d_factor
            self.length = self.length * gr_d_factor

        # Programmed cell death
        if (random.random() > 1.0 - mortility_rate):
            self.living = False

    def split(self):
        """
        split bacteria if critical length is reached, create new bacteria with new values
        :return: daughter bacteria
        """
        volumeRatio = 0.4 + 0.2 * random.random()
        volume = volumeRatio * self.getVolume()
        # Additional growth loss
        # ------looks better-------
        # if(True):
        # volume = volume*0.8*0.8*0.8

        # Advanced new position: random angular component,
        #                       radial component sum of two radii*0.8
        spawn_angle = [self.angle[0] + (0.5 - random.random()) * 3.1415 * 0.5,
                       self.angle[1] + (0.5 - random.random()) * 3.1415 * 0.5]
        self.angle = [self.angle[0] + (0.5 - random.random()) * 3.1415 * 0.5,
                      self.angle[1] + (0.5 - random.random()) * 3.1415 * 0.1]
        # Small probability: also rotate cell orientation by 180°
        if (random.random() > 0.8):
            spawn_angle[0] = spawn_angle[0] + 3.1415
        if (random.random() > 0.8):
            self.angle[0] = self.angle[0] + 3.1415
        # spawn_distance = (self.getVolume())**(1/3)*0.75
        spawn_distance = (self.length) * 1.0 * 0.1
        # print(self.getVolume())
        # print(spawn_distance)
        spawn_offsetX = spawn_distance * math.sin(spawn_angle[0]) * math.cos(spawn_angle[1])
        spawn_offsetY = spawn_distance * math.cos(spawn_angle[0]) * math.cos(spawn_angle[1])
        spawn_offsetZ = spawn_distance * math.sin(spawn_angle[1])
        # Move cell with the ratio 1/mass
        # massFactor = (self.getVolume()-volume)/(self.getVolume())
        pos = (self.pos[0] + spawn_offsetX * volumeRatio, self.pos[1] + spawn_offsetY * volumeRatio,
               self.pos[2] + spawn_offsetZ * volumeRatio)
        # print(str(self.pos) + " " + str((spawn_offsetX, massFactor, spawn_offsetY)))
        # bac = self(pos, (0,0,0), volume, 0)
        # width  = self.width#0.5*volume**(1/3)
        length = volumeRatio * self.length  # (volume*volumeRatio)**(2/3)/width
        bac = bacterium(pos, self.width, length, (0, 0, 0), spawn_angle, 0, True, False)
        # print(str(pos) +" "+ str(width) +" "+ str(length) +" "+ str((0,0,0)) +" "+ str(random.random()*3.14))

        # Also modify parent self
        # self.width  = self.width #*(1-volumeRatio)**(1/3)*0.8#growth_factor
        # self.length = self.length*(1-volumeRatio)**(1/3)#*growth_factor
        self.length = (1 - volumeRatio) * self.length  # (volume*(1-volumeRatio))**(2/3)/width

        self.pos[0] = self.pos[0] - spawn_offsetX * (1.0 - volumeRatio)
        self.pos[1] = self.pos[1] - spawn_offsetY * (1.0 - volumeRatio)
        self.pos[2] = self.pos[2] - spawn_offsetZ * (1.0 - volumeRatio)

        # return new bacterium
        return bac

    def getPositions(self):
        positions = []
        dx_length = self.length * math.sin(self.angle[0]) * math.cos(self.angle[1])
        dy_length = self.length * math.cos(self.angle[0]) * math.cos(self.angle[1])
        dz_length = self.length * math.sin(self.angle[1])
        positions.append((self.pos[0] - int(0.1 * dx_length), self.pos[1] - int(0.1 * dy_length),
                          self.pos[2] - int(0.1 * dz_length)))
        for index in range(-int(self.length * 0.01), int(self.length * 0.01) + 1):
            positions.append((self.pos[0] + int(10 * index * math.sin(self.angle[0]) * math.cos(self.angle[1])),
                              self.pos[1] + int(10 * index * math.cos(self.angle[0]) * math.cos(self.angle[1])),
                              self.pos[2] + int(10 * index * math.sin(self.angle[1]))))
        positions.append((self.pos[0] + int(0.1 * dx_length), self.pos[1] + int(0.1 * dy_length),
                          self.pos[2] + int(0.1 * dz_length)))
        return positions

    def getVolume(self):
        """ gives out the cubic volume equivalent """
        return self.width * self.width * self.length

    def move(self, frameDim, friction, dt):
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
        for position in self.getPositions():
            if (False):
                xborder = frameDim[1] * 0.5
                if (position[0] < -xborder + self.width):
                    self.velocity[0] = self.velocity[0] + (-position[0] - xborder + self.width) * 0.1
                    self.totalForce_equivalent = self.totalForce_equivalent + absolute(
                        -position[0] - xborder + self.width) * 0.1
                if (position[0] > xborder - self.width):
                    self.velocity[0] = self.velocity[0] + (-position[0] + xborder - self.width) * 0.1
                    self.totalForce_equivalent = self.totalForce_equivalent + absolute(
                        -position[0] + xborder - self.width) * 0.1
            yborder = frameDim[0] * 0.5
            if (position[1] < -yborder + self.width):
                self.velocity[1] = self.velocity[1] + (-position[1] - yborder + self.width) * 0.1
                self.totalForce_equivalent = self.totalForce_equivalent + absolute(
                    -position[1] - yborder + self.width) * 0.1
            #    if(position[1]>yborder-self.width):
            #        self.velocity[1] = self.velocity[1] + (-position[1]+yborder-self.width)*0.1    
            #        self.totalForce_equivalent = self.totalForce_equivalent + absolute(-position[1]+yborder-self.width)*0.1    
            # Bottom-boundary (z<=0)
            if (position[2] < -0 + self.width):
                self.velocity[2] = self.velocity[2] + (-position[2] + self.width) * 0.1
                # self.totalForce_equivalent = self.totalForce_equivalent + absolute(-position[2]+self.width)*0.1

                # Bounding-box torque
                positions = self.getPositions()
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

        self.pos[0] = self.pos[0] + self.velocity[0] * dt
        self.pos[1] = self.pos[1] + self.velocity[1] * dt
        self.pos[2] = self.pos[2] + self.velocity[2] * dt
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
        self.totalForce_equivalent = self.totalForce_equivalent * 0.9

    def interaction(self, _Bacterium):
        """ return interaction term with bacteria in local environment"""
        lenPos = len(self.getPositions())

        dx = _Bacterium.pos[0] - self.pos[0]
        dy = _Bacterium.pos[1] - self.pos[1]
        dz = _Bacterium.pos[2] - self.pos[2]
        dr = dx * dx + dy * dy + dz * dz
        interactionfactor = 0.5
        # If the bacterium is "far away"
        # Ignore all circles of this Bacteria to improve speed
        # Do the following operation with all other Bacteria
        far_away_factor = 8
        if (dr ** 0.5 < (self.width) * 1.05 * 5 * far_away_factor):
            positions = self.getPositions()
            for index in range(lenPos - 1):
                position = positions[index]

                dx = _Bacterium.pos[0] - position[0]  # self.pos[0]
                dy = _Bacterium.pos[1] - position[1]  # self.pos[1]
                dz = _Bacterium.pos[2] - position[2]  # self.pos[2]
                dr = dx * dx + dy * dy + dz * dz
                interactionfactor = 0.5
                # If the bacterium is "far away"
                # Ignore all circles of this one to improve speed
                far_away_factor = 4
                if (dr ** 0.5 < (self.width) * 1.05 * 5 * far_away_factor):

                    # Each not only with the center of _Bacterium, instead also every sphere of this
                    _positions = _Bacterium.getPositions()
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
                            if (dr ** 0.5 < (self.width) * 1.05 * 5):
                                repulsion_x = -dx * 40 / dr  # (1.5)
                                repulsion_y = -dy * 40 / dr  # (1.5)
                                repulsion_z = -dz * 40 / dr  # (1.5)

                            # New repulsion-function design
                            # if(dr**0.5<(Bacterium.getVolume())**(1/3)*1.25*1.7):
                            #    repulsion_x = -dx*40    /dr**(0.5)*repulsion_function(Bacterium.getVolume(), dr)
                            #    repulsion_y = -dy*40    /dr**(0.5)*repulsion_function(Bacterium.getVolume(), dr)

                            #
                            self.velocity[0] = self.velocity[0] + int(
                                repulsion_x / lenPos ** (1 / 2) * interactionfactor)
                            self.velocity[1] = self.velocity[1] + int(
                                repulsion_y / lenPos ** (1 / 2) * interactionfactor)
                            self.velocity[2] = self.velocity[2] + int(
                                repulsion_z / lenPos ** (1 / 2) * interactionfactor)
                            # add torque
                            t_radius = (index - lenPos * 0.5)
                            self.velocity_angular[0] = self.velocity_angular[0] + t_radius * math.cos(
                                self.angle[0]) * math.cos(
                                self.angle[1]) * repulsion_x / lenPos * 0.05 * interactionfactor
                            self.velocity_angular[0] = self.velocity_angular[0] - t_radius * math.sin(
                                self.angle[0]) * math.cos(
                                self.angle[1]) * repulsion_y / lenPos * 0.05 * interactionfactor
                            # Torque on second angle
                            self.velocity_angular[1] = self.velocity_angular[1] + t_radius * math.cos(self.angle[1]) * (
                                        repulsion_x ** 2 + repulsion_y ** 2) ** (
                                                                   1 / 2) / lenPos * 0.0125 * interactionfactor
                            self.velocity_angular[1] = self.velocity_angular[1] + t_radius * math.sin(
                                self.angle[1]) * repulsion_z / lenPos * 0.05 * interactionfactor
                            self.totalForce_equivalent = self.totalForce_equivalent + absolute(
                                repulsion_x / lenPos * interactionfactor)
                            self.totalForce_equivalent = self.totalForce_equivalent + absolute(
                                repulsion_y / lenPos * interactionfactor)

                            # Actio-Reactio
                            _Bacterium.velocity[0] = _Bacterium.velocity[0] - (
                                        repulsion_x / lenPos ** (1 / 2) * interactionfactor)
                            _Bacterium.velocity[1] = _Bacterium.velocity[1] - (
                                        repulsion_y / lenPos ** (1 / 2) * interactionfactor)
                            _Bacterium.velocity[2] = _Bacterium.velocity[2] - (
                                        repulsion_z / lenPos ** (1 / 2) * interactionfactor)
                            _Bacterium.totalForce_equivalent = _Bacterium.totalForce_equivalent + absolute(
                                repulsion_y / lenPos * interactionfactor)
                            _Bacterium.totalForce_equivalent = _Bacterium.totalForce_equivalent + absolute(
                                repulsion_x / lenPos * interactionfactor)

        return _Bacterium, self

    # Actual not used structural part due to the interaction part
    # and resulting design difficulties
    # see core function
    '''
    def evolution(self, duration):
        """ evolution of cell for duration, break evolution if new cell is created  """
        iterations = []
        bacteria = [self]
        for t in range(0, duration):
            for bac in bacteria:
                bac.move()
                bac.grow()
                #if bac.length >= bac.split_length:
                #    bac, new_bac = bac.split()
                #    bacteria.extend((bac, new_bac))

            iterations.append(bacteria)

        return iterations
    '''

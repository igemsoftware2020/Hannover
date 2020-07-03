import numpy as np
from matplotlib import pyplot as plt
import random


class Bacteria(object):
    growth_rate = 0.2  # um / s , use Monod equation if differenz growth rates with substrate concentration
    split_length = 70   # cell splits, if height is greater that 10 um

    def __init__(self, coord, width, length, velocity, angle, index, **kwargs):
        """
        initialize a instance of the Bacteria class
        :param coord: position of bacteria center [x_pox, y_pos, z_pos]
        :param width: diameter of ellipse
        :param length:  length of ellipse
        :param velocity: velocity of bacteria [v_x, v_y, v_z]
        :param angle: angle of bacteria in radian  measured to x axis
        :param index: integer value, index of bacteria
        """

        self.coord = coord
        self.width = width    # bacteria is a ellipse with radii width / length
        self.length = length  # lenghts in um
        self.angle = angle    # angle is orientation of cell measured to x- axes in degree
        self.v_x, self.v_y = velocity[0], velocity[1]
        self.index = index

    def print_info(self):
        """
        print all values of parameters of bacteria
        :return:
        """
        print("\n**** Current values of bacci {index} ****\n"
              "  Position       ({x_pos:f}, {y_pos:f})\n"
              "  Width, Length      {width:.2f},{length:.2f}\n"
              "  angle            {o:.1f}\n"
              "  velocity           {vx:.2f}, {vy:.2f}".format(index=self.index, x_pos=self.coord[0],
                                                        y_pos=self.coord[1],
                                                        width=self.width, length=self.length,
                                                        o=self.angle, vx=self.v_x, vy=self.v_y)
              )

    def grow(self):
        """
        grow bacteria for 1 second with speed growth_rate
        """
        self.width = self.width * (1 + self.growth_rate / 5)
        self.length = self.length * (1 + self.growth_rate)

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

            angle_to_y = np.arccos(np.dot(center, np.array([0, 1])) / (np.linalg.norm(center))) # top case
            if center[0] < 0:
                angle_to_y = -angle_to_y
            c1 = center
            c2 = center + distance * rotate_vector(rotate_vector(center, angle), angle_to_y)
            return c1, c2
        if self.length >= self.split_length:
            # update parameters of mother
            self.length = self.length / 2
            self.width = self.width / 2
            self.v_x, self.v_y = self.v_x / 2, self.v_y / 2
            # calculate new centers
            c1, c2 = new_centers(self.coord, self.angle, 2*self.length)
            # set new bacteria
            self.coord = c1
            bac_1 = Bacteria(coord=c1, width=self.width, length=self.length, velocity=[self.v_x, self.v_y],
                            angle=self.angle, index=self.index)
            bac_2 = Bacteria(coord=c2, width=self.width, length=self.length, velocity=[-self.v_x, -self.v_y],
                            angle=self.angle, index=self.index + 1)
            del self
            return bac_1, bac_2

    def move(self):
        """
        Move Bacteria for 1 time unit and add random movement and rotation
        :return:
        """
        rand_x, rand_y = brownian(100)
        self.coord = np.array([self.coord[0] + self.v_x * 1 + rand_x[0],
                               self.coord[1] + self.v_y * 1 + rand_y[0]])

        self.angle = self.angle * np.random.normal(loc=1, scale=2, size=1)[0]

    def interaction(self):
        """ return interaction term with bacteria in local environment"""

    def evolution(self, duration):
        """ evolution of cell for duration, break evolution if new cell is created  """
        iterations = []
        bacteria = [self]
        for t in range(0, duration):
            for bac in bacteria:
                bac.move()
                bac.grow()
                if bac.length >= bac.split_length:
                    bac, new_bac = bac.split()
                    bacteria.extend((bac, new_bac))

            iterations.append(bacteria)

        return iterations


def plot_ellipse(x_center, y_center, width, height, angle):
    """ function returning x and y values of an ellipse for plotting"""
    t = np.linspace(0, 2*np.pi, 1000)
    # ellipse = np.array([x_center + width * np.cos(t), y_center + height * np.sin(t)])
    ellipse_rotated = np.array([x_center + width * np.cos(t) * np.cos(angle) - height * np.sin(t) * np.sin(angle),
                                y_center + width * np.cos(t) * np.sin(angle) + height * np.sin(t) * np.cos(angle)])
    return ellipse_rotated[0], ellipse_rotated[1]


def plot_bacteria(bacteria):
    """ plots all bacteria in list bacteria as ellipses"""
    bac_plots = []
    bac_centers = []
    for bac in bacteria:
        x, y = plot_ellipse(x_center=bac.coord[0], y_center=bac.coord[1],
                            width=bac.width, height=bac.length, angle=bac.angle)
        x_pos, y_pos = bac.coord[0], bac.coord[1]
        bac_plots.append([x, y])
        bac_centers.append([x_pos, y_pos])

    ax = plt.subplot()
    for bac_plot, bac_center in zip(bac_plots, bac_centers):
        ax.plot(bac_plot[0], bac_plot[1], color="grey")
        ax.scatter(bac_center[0], bac_center[1])

    plt.show()


def brownian(n):
    """ simulating random movement for n steps, returns 2 times n random numbers"""
    x = np.cumsum(np.random.randn(n))
    y = np.cumsum(np.random.randn(n))
    return x, y


def test_evolution():
    initial_vals = {"position": np.array([4, 4]), "width": 15, "length": 30,
                    "velocity": np.array([10, -10]), "angle": np.pi / 4, "index": 0}
    bac = Bacteria(coord=initial_vals["position"], width=initial_vals["width"],
                   length=initial_vals["length"], velocity=initial_vals["velocity"],
                   angle=initial_vals["angle"], index=initial_vals["index"])
    iterations = bac.evolution(9)


    frame = [np.ones((200, 200, 3)) for i in iterations]
    for i in range(0, len(iterations)):
        for bac in iterations[i]:
            cv2.ellipse(image=frame[i], centerCoordinates=bac.coord, axesLength=(bac.length, bac.width), angle=bac.angle, thickness=2)

    plt.imshow(frame[1])


def test_splitting_moving():
    initial_vals = {"position": np.array([10, 10]), "width": 2, "length": 5,
                    "velocity": np.array([0.1, 0.05]), "angle": np.pi / 4, "index": 0}
    bac = Bacteria(coord=initial_vals["position"], width=initial_vals["width"], length=initial_vals["length"], velocity=initial_vals["velocity"], angle=initial_vals["angle"], index=initial_vals["index"])

    while bac.length < bac.split_length:
        bac.move()
        bac.grow()

    bac, new_bac = bac.split()

    x_val, y_val = plot_ellipse(x_center=bac.coord[0], y_center=bac.coord[1], width=bac.width, height=bac.length, angle=bac.angle)
    x_val_1, y_val_1 = plot_ellipse(new_bac.coord[0], new_bac.coord[1], new_bac.width, new_bac.length, new_bac.angle)

    ax = plt.subplot()
    ax.plot(x_val, y_val, 'blue')
    ax.plot(x_val_1, y_val_1, "orange")
    ax.scatter(new_bac.coord[0], new_bac.coord[1], alpha=0.5,c='orange' )
    ax.scatter(bac.coord[0], bac.coord[1], c='blue')
    ax.grid(color='lightgray', linestyle='--')
    ax.set_xlim([0,30])
    ax.set_ylim([0,30])
    ax.axis("equal")
    plt.show()


test_evolution()


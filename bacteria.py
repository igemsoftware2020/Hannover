import numpy as np
from matplotlib import pyplot as plt


def brownian_movement():
    np.random.seed(1234)
    """
    Simulates a Brownian motion
    :param int N : the number of discrete steps
    :param int T: the number of continuous time steps
    :param float h: the variance of the increments
    """
    def random_movement(N, T, h):
        dt = 1. * T / N  # total number of time steps
        random_increments = np.random.normal(0.0, 1.0 * h, N) * np.sqrt(dt)  # the epsilon values
        brownian_motion = np.cumsum(random_increments)  # calculate the brownian motion
        # brownian_motion = np.insert(brownian_motion, 0, 0.0)  # insert the initial condition
        return brownian_motion, random_increments

    N = 50  # the number of discrete steps
    T = 1  # the number of continuous time steps
    h = 5  # the variance of the increments
    # generate a brownian motion
    X, epsilon = random_movement(N, T, h)
    return X


class Bacteria(object):
    growth_rate = 0.01  # um / s , use Monod equation if differenz growth rates with substrate concentration
    split_length = 15   # cell splits, if height is greater that 10 um

    def __init__(self, coord, width, length, velocity, orientation, index):
        self.coord = coord
        self.width = width  # bacteria is a ellipse with radii width / height
        self.length = length  # lenghts in um
        self.orientation = orientation  # orientation is angle measured to x- axes in degree
        self.v_x, self.v_y = velocity[0], velocity[1]
        self.index = index

    def print_info(self):
        print("\n**** Current values of bacci {index} ****\n"
              "  Position       ({x_pos:f}, {y_pos:f})\n"
              "  Width, Length      {width:.2f},{length:.2f}\n"
              "  Orientation            {o:.1f}\n"
              "  velocity           {vx:.2f}, {vy:.2f}".format(index=self.index, x_pos=self.coord[0],
                                                        y_pos=self.coord[1],
                                                        width=self.width, length=self.length,
                                                        o=self.orientation, vx=self.v_x, vy=self.v_y)
              )

    def grow(self):
        self.width = self.width * (1 + self.growth_rate / 6)
        self.length = self.length * (1 + self.growth_rate)

    def split(self):
        print("***** Splitting *****")
        if self.length >= self.split_length:
            self.length = self.length / 2
            self.width = self.width / 2
            self.v_x = self.v_x / 2; self.v_y = self.v_y / 2
            x, y = self.coord
            x1, y1 = x + self.length * np.cos(self.orientation*180/np.pi), \
                     y + self.length * np.sin(self.orientation*180/np.pi)
            x2, y2 = x - self.length * np.cos(self.orientation*180/np.pi),\
                     y - self.length * np.sin(self.orientation*180/np.pi)
            self.coord = [x2, y2]

            new_bac = Bacteria(coord=[x1, y1], width=self.width,
                               length=self.length, velocity=[self.v_x, self.v_y],
                               orientation=self.orientation, index=self.index+1)
            return new_bac

    def move(self):
        self.coord = self.coord[0] + self.v_x * 1, self.coord[1] + self.v_y * 1

    def evolution(self):
        """ TODO"""
        pass


bac = Bacteria(coord=[0, 0], width=1, length=3, velocity=[0, 0], orientation=0, index=0)
pos_x = []
pos_y = []
while bac.length < bac.split_length:
    bac.grow()
    bac.move()
    pos_x.append(bac.coord[0]); pos_y.append(bac.coord[1])

new_bac = bac.split()

plt.scatter(new_bac.coord[0], new_bac.coord[1])
plt.scatter(pos_y, pos_y)
plt.show()

brownian_movement()

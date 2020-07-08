import numpy as np
import matplotlib.pyplot as plt
import bacteria


class Biofilm(object):

    def __init__(self):
        self.living_bacteria = []
        self.evolution = []

    def initial_configuration(self, initial):
        """
        Set initial configuration of Biofilm

        """
        if initial is None:
            initial = {"position": np.array([4, 4]), "length": 30,
                       "velocity": np.array([10, -10]), "angle": np.pi / 4, "index": 0}

        bac = bacteria.Bacteria(pos=initial["position"], length=initial["length"],
                                velocity=initial["velocity"], angle=initial["angle"],
                                index=initial["index"])
        self.living_bacteria = [bac]
        self.evolution.append([bac])

    def print_current_state(self):
        print("***** Current state of Biofilm ******\n"
              "     Bacteria alive : {alive:.0f}\n"
              "     Evolution step No. {step:.0f}".format(alive=len(self.living_bacteria),
                                                          step=len(self.evolution))
              )

    def evolve_bacteria(self):
        evolved = []
        for bac in self.living_bacteria:
            bac.move()
            bac.grow()
            if bac.length >= bac.split_length:
                bac, new_bac = bac.split()
                evolved.extend((bac, new_bac))
            else:
                evolved.append(bac)
        self.evolution.append(evolved)
        self.living_bacteria = evolved

    def plot_iterations(self):
        i = 0
        for steps in self.evolution:
            print("Currently depicting step {} ".format(i))
            bac_plots = []
            bac_centers = []
            for bac in steps:
                x, y = Biofilm.plot_ellipse(x_center=bac.pos[0], y_center=bac.pos[1],
                                    width=bac.width, height=bac.length, angle=bac.angle)
                x_pos, y_pos = bac.pos[0], bac.pos[1]
                bac_plots.append([x, y])
                bac_centers.append([x_pos, y_pos])

            ax = plt.subplot()
            for bac_plot, bac_center in zip(bac_plots, bac_centers):
                ax.plot(bac_plot[0], bac_plot[1], color="grey")
                ax.scatter(bac_center[0], bac_center[1])
            plt.show()
            i += 1

    @staticmethod
    def plot_ellipse(x_center, y_center, width, height, angle):
        """ function returning x and y values of an ellipse for plotting"""
        t = np.linspace(0, 2*np.pi, 1000)
        ellipse_rotated = np.array([x_center + width * np.cos(t) * np.cos(angle) - height * np.sin(t) * np.sin(angle),
                                   y_center + width * np.cos(t) * np.sin(angle) + height * np.sin(t) * np.cos(angle)])
        return ellipse_rotated[0], ellipse_rotated[1]



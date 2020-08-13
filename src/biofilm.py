# coding latin-1


# ********************************************************************************************
# imports
import numpy as np
import random

# custom libraries
from src.bacteria import Bacterium


class Biofilm(object):

    def __init__(self):
        self.bacteria = []

    def spawn(self, startNum):
        for _ in range(startNum):
            pos = ((random.random() - 0.5) * 400.0, (random.random() - 0.5) * 400.0, 0)
            # volume = 3300+1700*random.random()
            width = 0.4 * 15  # *volume**(1/3)
            # width  = 0.4 *15#*volume**(1/3)
            length = 4 * 8 * 15  # *volume**(1/3)
            bac = Bacterium(pos, width, length, (0, 0, 0), [random.random() * 3.14, (0.5 - random.random() * 3.14)], 0,
                            True, False)
            self.bacteria.append(bac)

    def print_current_state(self):
        print("***** Current state of Biofilm ******\n"
              "     Bacteria alive : {alive:.0f}\n"
              "     ".format(bacteria=len(self.bacteria))
              )

    def evolve(self, gr_factor, gr_pr_i, gr_factor_inv, gr_d_factor, mortality_rate, motion_activation,
               motion_deactivation, frameDim, combined_friction, splitting_rate, dt):
        for Bacterium in self.bacteria:
            Bacterium.grow(gr_factor, gr_pr_i, gr_factor_inv, gr_d_factor, mortality_rate)

            if (Bacterium.living == False) and (Bacterium.getVolume() < 3300):
                self.bacteria.remove(Bacterium)
            # Manage repulsion
            # (The most expensive task)
            for _Bacterium in self.bacteria:
                if ((Bacterium.pos[0] == _Bacterium.pos[0]) and (Bacterium.pos[1] == _Bacterium.pos[1]) and (
                        Bacterium.pos[2] == _Bacterium.pos[2])):
                    True
                else:
                    [_Bacterium, Bacterium] = Bacterium.interaction(_Bacterium)

            Bacterium.move(frameDim, combined_friction, dt)

            # Manage Bacterial splitting
            # Add a little bit of random until -----it looks good and real-----
            if (Bacterium.getVolume() > 20000) and (Bacterium.living) and (Bacterium.moving == False):
                if (random.random() > splitting_rate) or (Bacterium.getVolume() > 35000):
                    # Split Bacterium into two,
                    # Differate between old and new cell
                    # Now add new bacterium to
                    self.bacteria.append(Bacterium.split())
            # Manage Bacterial Motion-Mode
            # Enable/Disable motion mode
            if (Bacterium.getVolume() > 20000) and (Bacterium.living):
                if (random.random() > 1.0 - motion_activation):
                    Bacterium.moving = True
                if (random.random() > 1.0 - motion_deactivation):
                    Bacterium.moving = False
            else:
                Bacterium.moving = False

    # This function should be removed here and included in core library
    # If using OpenGL with depth sorting enabled,
    # This function also gets redundant
    def sortbydepth(self, axis, _reverse):
        sorted_bacteria = self.bacteria
        # To return a new list, use the sorted() built-in function...
        return sorted(sorted_bacteria, key=lambda x: x.pos[axis], reverse=_reverse)

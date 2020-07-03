#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a awesome
        python script!
"""

import numpy as np
import matplotlib.pyplot as plt
import cv2

# own modules
import bacteria


def test_evolution():
    initial_vals = {"position": np.array([4, 4]), "width": 15, "length": 30,
                    "velocity": np.array([10, -10]), "angle": np.pi / 4, "index": 0}
    bac = bacteria.Bacteria(coord=initial_vals["position"], width=initial_vals["width"],
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
    bac = bacteria.Bacteria(coord=initial_vals["position"], width=initial_vals["width"], length=initial_vals["length"], velocity=initial_vals["velocity"], angle=initial_vals["angle"], index=initial_vals["index"])

    while bac.length < bac.split_length:
        bac.move()
        bac.grow()

    bac, new_bac = bac.split()

    x_val, y_val = bacteria.plot_ellipse(x_center=bac.coord[0], y_center=bac.coord[1], width=bac.width, height=bac.length, angle=bac.angle)
    x_val_1, y_val_1 = bacteria.plot_ellipse(new_bac.coord[0], new_bac.coord[1], new_bac.width, new_bac.length, new_bac.angle)

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

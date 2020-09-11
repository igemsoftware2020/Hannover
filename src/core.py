# coding latin-1
# ********************************************************************************************
# THR
# semiautodeveloped code
# ********************************************************************************************


# ********************************************************************************************
# imports

import math
import os

#mport cv2
import numpy as np

from biofilm import Biofilm
from constants import Constants as C
from utils import plot_size, plot_force, plot_velocities, plot_positions, bacteria_as_pandas, get_info_file_path, \
    prompt_log_at_start, plot_num


# custom libraries


# ********************************************************************************************
# general dependency functions
def norm(value):
    if value <= 1.0:
        if value >= 0.0:
            return value
        else:
            return -value
    else:
        return 1.0


# ********************************************************************************************
# core function


# ********************************************************************************************
# Configuration of the biofilm
frameDim = (500, 750)
startNum = 1
# Medium specification
combined_friction = 0.63

# living cells specifications
gr_factor = 1.1  # 1.0
gr_factor_inv = 0.05  # 0.2
gr_pr_i = 10.1  # 1.1
motion_activation = 0.005
motion_deactivation = 0.01
mortility_rate = 0.0006
splitting_rate = 0.96

# not living cells specifications
gr_d_factor = 0.999


def video_init():
    # get path to python code to retrieve parent folder
    # save to output directory
    path = os.path.dirname(os.path.realpath(__file__))
    # parentDir = "\\".join(path.split("\\")[:-1])
    ROOT_PATH = os.getcwd()
    # outputDir = parentDir + "\\output\\"
    outputDir = os.path.join(ROOT_PATH, "output\\")
    print("outputDir" + outputDir)
    # Add filename
    num = 1
    path_out = outputDir + "core_test1_1_" + str(num) + ".avi"
    while (os.path.exists(path_out)):
        num = num + 1
        path_out = outputDir + "core_test1_1_" + str(num) + ".avi"

    fourcc = cv2.VideoWriter_fourcc(*'XVID')
    out = cv2.VideoWriter(path_out,
                          fourcc,
                          60.0, (frameDim[1], frameDim[0]))
    return (out, path_out)


def storage_file_init(video_path):
    # use corresponding video file path
    # by changing .avi to .txt
    path_out = video_path.replace(".avi", ".txt")
    # create simple text file
    txt = open(path_out, "w")
    # Add info line to text file
    txt.write("#********************************************************************************************\n")
    txt.write("#THR\n")
    txt.write("#This file contains all relevant parameter of the corresponding video\n")
    txt.write("#********************************************************************************************\n\n")

    txt.write("#coding:\n")
    txt.write("#index, totalForce_equivalent, living, moving, width, length,\n")
    txt.write("#[pos[0], pos[1], pos[2]], [velocity[0], velocity[1], velocity[2]],\n")
    txt.write("#[angle[0], angle[1]], [velocity_angular[0], velocity_angular[1]];\n")
    txt.write("#********************************************************************************************\n")
    return txt


def storage_file_append_frame(txt, num_frame):
    txt.write("\nframe " + str(num_frame) + "\n")


def storage_file_append_bacterium(txt, bacterium, index):
    # txt.write("#index; totalForce_equivalent; living; moving; width; length;\n")
    txt.write(str(index) + "; " + str(bacterium.total_force) + "; " + str(bacterium.living) + "; ")
    txt.write(str(bacterium.moving) + "; " + str(bacterium.width) + "; " + str(bacterium.length) + "; ")
    # txt.write("#[pos[0]; pos[1]; pos[2]]; [velocity[0]; velocity[1]; velocity[2]];\n")
    txt.write(str([bacterium.position[0], bacterium.position[1], bacterium.position[2]]) + "; ")
    txt.write(str([bacterium.velocity[0], bacterium.velocity[1], bacterium.velocity[2]]) + "; ")
    # txt.write("#[angle[0]; angle[1]]; [velocity_angular[0]; velocity_angular[1]];\n")
    txt.write(str([bacterium.angle[0], bacterium.angle[1]]) + "; ")
    txt.write(str([bacterium.velocity_angular[0], bacterium.velocity_angular[1]]) + ";\n")


# ********************************************************************************************
def coreFunction():
    # Initialization
    (out, path_out) = video_init()
    (txt) = storage_file_init(path_out)

    # Seed
    biofilm = Biofilm()
    biofilm.spawn()
    # coreLoop
    coreLoop(biofilm, out, txt)


def coreLoop(biofilm, out, txt):
    # The frames are only counted for the data storage file
    num_frame = 0
    storage_file_append_frame(txt, num_frame)
    # Loop until the frame is closed
    cv2.imshow("frame", 5)
    # This unbeatiful call of the frame has to be done
    # In order to exit on exit button the way below
    while (cv2.getWindowProperty("frame", 0) >= 0):

        # Timestep
        dt = 1
        # (should be directly taken from system to improve visualization)
        #   ----correction: since this is calculated and prestored in video file with fixed framerate,
        #                   a constant preset timedifference is optimal
        frameDimBacteria = (int(frameDim[0]), int(frameDim[1]))  # Test condition field_size=*0.25
        biofilm.evolve()

        # Draw the frame
        # Dark_simmulation
        frame = np.zeros((frameDim[0], frameDim[1], 3))
        frame[:] = (0.0, 0.2, 0.2)

        # Sort bacteria by depth
        # Has to be done right now for correct depth visualization
        bacteria = biofilm.sort_by_depth(2, False)
        bacterium_index = 0
        for Bacterium in bacteria:
            bacterium_index = bacterium_index + 1
            # Color coding of the overall-force-sum
            positions = Bacterium.get_position()
            for pos_index in range(len(positions)):
                if (pos_index < len(positions)):
                    pos = positions[pos_index]
                else:
                    pos = positions[0]  # Draw the zero_circle at last for correct bounding lines
                # Bacteria_height coded by cell-brightness
                z_coding = norm(100.0 / (100.0 + 50 - pos[2] + Bacterium.width + 0.001))
                # *************************************************************************
                # Front-View
                pos_0 = 0 + pos[0]
                pos_1 = 0 + pos[1]
                _angle = Bacterium.angle[0] * 180 / 3.14

                rotated = True
                if (math.cos(Bacterium.angle[1]) < 0):
                    rotated = False

                if (Bacterium.living == True):
                    # b = norm(255 * (1.0 - 1.0 * 750.0 / (750.0 + Bacterium.total_force)))
                    b = 1E-5
                    # frame = cv2.circle(frame, (int(frameDim[1] * 0.5 + pos_0), int(frameDim[0] * 0.5 + pos_1)),
                    # int((Bacterium.width * 2)),
                    # (1.0, norm((0.4 + b * 0.4 + 0.1 * z_coding)), norm(0.2 + 0.2 * z_coding)), -1)
                else:
                    # b = 1 * (1.0 - 1.0 * 250.0 / (250.0 + Bacterium.total_force))
                    b = 1E-5
                    # frame = cv2.circle(frame, (int(frameDim[1] * 0.5 + pos_0), int(frameDim[0] * 0.5 + pos_1)),
                    # int((Bacterium.width * 2)), (0.3, (0.25 + b), (0.25 + b)), -1)

                contour_color = (0.5, 0.4, 0.2)
                frame = cv2.ellipse(frame, (int(frameDim[1] * 0.5 + pos_0), int(frameDim[0] * 0.5 + pos_1)),
                                    (int((Bacterium.width * 2)), int((Bacterium.width * 2))),
                                    0, -20 - _angle, 20 - _angle, contour_color, 2)
                frame = cv2.ellipse(frame, (int(frameDim[1] * 0.5 + pos_0), int(frameDim[0] * 0.5 + pos_1)),
                                    (int((Bacterium.width * 2)), int((Bacterium.width * 2))),
                                    0, 180 - 20 - _angle, 180 + 20 - _angle, contour_color, 2)
                if pos_index == len(positions) - 1:
                    if rotated == False:
                        frame = cv2.ellipse(frame, (int(frameDim[1] * 0.5 + pos_0), int(frameDim[0] * 0.5 + pos_1)),
                                            (int((Bacterium.width * 2)), int((Bacterium.width * 2))),
                                            0, 20 - _angle, -180 + 20 - _angle, contour_color, 2)
                    else:
                        frame = cv2.ellipse(frame, (int(frameDim[1] * 0.5 + pos_0), int(frameDim[0] * 0.5 + pos_1)),
                                            (int((Bacterium.width * 2)), int((Bacterium.width * 2))),
                                            0, 180 - 20 - _angle, 20 - _angle, contour_color, 2)

                if pos_index == 0:
                    if rotated == False:
                        frame = cv2.ellipse(frame, (int(frameDim[1] * 0.5 + pos_0), int(frameDim[0] * 0.5 + pos_1)),
                                            (int((Bacterium.width * 2)), int((Bacterium.width * 2))),
                                            0, 180 - 20 - _angle, 20 - _angle, contour_color, 2)
                    else:
                        frame = cv2.ellipse(frame, (int(frameDim[1] * 0.5 + pos_0), int(frameDim[0] * 0.5 + pos_1)),
                                            (int((Bacterium.width * 2)), int((Bacterium.width * 2))),
                                            0, 20 - _angle, -180 + 20 - _angle, contour_color, 2)

            storage_file_append_bacterium(txt, Bacterium, bacterium_index)

        cv2.imshow("frame", frame)

        frame = cv2.normalize(frame, None, alpha=0, beta=255, norm_type=cv2.NORM_MINMAX, dtype=cv2.CV_32F)
        frame = frame.astype(np.uint8)

        # Store frame and data in corresponding files
        num_frame = num_frame + 1
        out.write(frame)
        if (bacterium_index > 0):
            storage_file_append_frame(txt, num_frame)

        key = cv2.waitKey(5)
        if (key != -1):
            if (key == 32):
                continue
            print(key)
            break
            out.release()
            txt.close()


def blind_run():
    info_file_name = get_info_file_path()
    info_file_path = info_file_name.parent

    print(prompt_log_at_start(info_file_name))

    biofilm = Biofilm()
    biofilm.simulate(duration_in_min=60, save_name=info_file_name)

    data = bacteria_as_pandas(info_file_name)

    plot_velocities(data, info_file_path, save_fig=True)
    plot_positions(data, info_file_path, save_fig=True)
    plot_force(data, info_file_path, save_fig=True)
    plot_size(data, info_file_path, save_fig=True)
    plot_num(data,info_file_path,save_fig=True )
    

def plot_testing(info_file_name):
    info_file_path = C.OUTPUT_PATH / info_file_name

    data = bacteria_as_pandas(info_file_path)
    plot_velocities(data, info_file_path)
    plot_positions(data, info_file_path)
    plot_force(data, info_file_path)
    plot_size(data, info_file_path)


# ********************************************************************************************
# main-method to start the program
# ********************************************************************************************
if __name__ == "__main__":
    blind_run()
    # coreFunction())

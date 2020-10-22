#!/usr/bin/env python
# -*- coding: utf-8 -*-
# imports
import json
import os
import tkinter as tk
from pathlib import Path
from tkinter import filedialog
from typing import Dict

import numpy as np
import pandas as pd
# custom libraries
from BiofilmSimulation.constants import Constants
from BiofilmSimulation.formulas import get_euclid_norm


# ********************************************************************************************
# Data handling functions for storing the data as json and reading from the json file



def write_log_template(info_file_path, constants: Constants):
    """ saves a json template for saving bacteria parameters"""
    with open(info_file_path, 'w+') as json_file:
        data = {'BACTERIA': {}, 'CONSTANTS': {}}
        data['CONSTANTS'].update(constants.get_simulation_constants())
        json.dump(data, json_file, indent=4)


def read_in_log(info_file_path) -> Dict:
    """ returns data from a json file as dictionary"""
    with open(info_file_path, "r") as json_file:
        data = json.load(json_file)
    return data


def ask_for_log_dir():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(initialdir=os.getcwd())
    return Path(file_path)


def bacteria_as_pandas(info_file_path) -> pd.DataFrame:
    """ sorts bacteria parameters into DataFrame"""
    data = read_in_log(info_file_path)
    return pd.DataFrame(data['BACTERIA']).transpose()


def constants_as_pandas(info_file_path):
    """ sorts all saved constant into a Dataframe and returns"""
    data = read_in_log(info_file_path)
    return pd.DataFrame(data['CONSTANTS']).transpose()


def save_dict_as_json(data: Dict, info_file_path: Path):
    with open(str(info_file_path), 'w') as json_file:
        json.dump(data, json_file)


def get_data_to_parameter(data: pd.DataFrame, key: str, exact: bool = False):
    """
    Gets complete bacteria data as a DataFrame.
    Resorts data to parameters 'key' into another DataFrame:
        bacteria_0_parameter    ... bacteria_N_parameter
    0       100                         NaN
    1       12                          NaN
    2       13                          10
    .
    .
    .

    If the key is position or velocity, calculates the absolute value first.
    Returned Data is sorted according to the iteration step.
    """
    dic = {}
    for index, bacteria in data.iterrows():
        df_bac = pd.DataFrame(bacteria).loc[key, :]
        for vectors in df_bac:
            if (key == 'velocity' or key == 'position' or key == 'acceleration') & (exact is False):
                # calculate norm for velocity and position
                vectors = get_euclid_norm(vectors)
            dic.update({str(index) + '_' + key: vectors})
    df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in dic.items()]))

    def isnan(vector):
        if np.isnan(np.min(vector)) and exact:
            return False
        elif exact:
            return True
        elif not exact:
            return pd.isna(vector)

    df = df.transform(lambda x: sorted(x, key=isnan, reverse=True))
    return df


def get_z(data):
    x = []
    y = []
    z = []
    for bac in data['position'].index:
        x.append(data['position'][bac][-1][0])
        y.append(data['position'][bac][-1][1])
        z.append(data['position'][bac][-1][2])
    a = np.vstack((x, y, z)).T
    pos = pd.DataFrame(a, columns=['x', 'y', 'z'])
    pos = pos.sort_values('z', ascending=False)
    pos = pos[0:int(np.round((len(pos) * 0.1)))].values
    x = pos[:, 0]
    y = pos[:, 1]
    z = pos[:, 2]
    return x, y, z

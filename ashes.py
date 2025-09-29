'''
See BEM-comparison-outputs for plots results of comparison

This file opens the raw data from the Ashes simulations of loads on the blades.
In the Ashes sensor panel, 'Thrust force, distr.' and 'Torque force distr.' are the ones
I've extracted. The raw output has a lot of crap from all the other sensors.
I wrote remove_header() function to shave off the upper header rows of the .txt output files.
The Extracted outputs from Ashes is stored in the 2D arrays p_n and p_t, starting with 5 m/s, then 9 and so forth.
The corresponding radii on the blade are found in vector 'span'.

The simulations in Ashes have been done at 0 degree pitch, Tip Speed Ratio = 8, and [5, 9, 11, 20] m/s
'''

import pandas as pd
import matplotlib.pyplot as plt
import os
from windtools import Tools

class AshesData:
    def __init__(self):
        self.windspeeds = [5, 9, 11, 20] # m/s
        self.p_n = [0] * 4 # [[p_n 5m/s], [p_n 9 m/s], ...]
        self.p_t = [0] * 4 # same

        self.BEM_p_n = []
        self.BEM_p_t = []
        self.BEM_rs = []

        #Blade span
        self.span = [0, 2.64302, 5.37977, 8.20274, 11.1031, 14.071, 17.0953, 20.1641, 23.2647, 26.3837, 29.5076, 32.6228, 35.7156, 38.773, 41.7824, 44.732, 47.6111, 50.4099, 53.1201, 55.7344, 58.247, 60.6534, 62.9501, 65.1352, 67.2076, 69.1675, 71.0159, 72.7545, 74.386, 75.9133, 77.3402, 78.6705, 79.9085, 81.0585, 82.1252, 83.113, 84.0265, 84.8703, 85.6487, 86.366]

    def remove_header(self, input_file, output_file):
        n_lines = 23;

        with open(input_file, "r") as f_in:
            lines = f_in.readlines()

        # Skip the first n_lines
        remaining = lines[n_lines:]

        with open(output_file, "w") as f_out:
            f_out.writelines(remaining)

    def compute_BEM(self, V0, TSR):
        lam = TSR # tip speed ratio
        theta_p = 0 # pitch angle
        R = 89.166  # m
        omega = lam * V0 / R
        rs, pts = [], []
        for row in range(1, len(pd.read_csv("bladedat.txt", sep=r"\s+", header=None)) - 1):
            r, c, beta, tc, R, fpath = Tools.read_blade_data(row)
            pn, pt, a, aprime, F = Tools.bem_single_element(r, c, beta, V0, omega, theta_p, fpath, R)
            rs.append(float(r)); pts.append(float(pt))

        self.BEM_p_n = pn
        self.BEM_p_t = pts
        self.BEM_rs = rs

    def extract_loads(self):
        files = ["5ms0deg8tps.txt",
                 "9ms0deg8tps.txt",
                 "11ms0deg8tps.txt",
                 "20ms0deg8tps.txt"]
        for file in files:
            file_name = os.getcwd() + "\\ashes-raw\\" + file
            out = os.getcwd() + "\\ashes-headers-removed\\" + file
            self.remove_header(file_name, out)

        # Read the file, skipping non-numeric header lines
        for i, file in enumerate(files):
            file_name = os.getcwd() + "\\ashes-headers-removed\\" + file
            df = pd.read_csv(
                file_name,
                sep=r"\s+",   # handles tabs and multiple spaces
                comment='#',             # if some lines start with '#'
                header=None,             # no header row
                engine="python"
            )
            #thrust force, torque force (distributed)
            end = df.shape[1]
            n_rows = df.shape[0]
            self.p_n[i] = df.iloc[n_rows - 1, end-(2*40):end - 40 + 1].values.tolist()
            self.p_t[i] = df.iloc[n_rows - 1, end-(40):end + 1].values.tolist()

# Program ->
ash = AshesData()       # instance of AshesData class
ash.extract_loads()     # clean up raw Ashes output, store distributed load data

for i in range(4):
    V0 = ash.windspeeds[i]      # Loop over wind speeds
    TSR = 8                     # TSR = 8 used in Ashes simulation
    ash.compute_BEM(V0, TSR)    # Compute BEM of blade with same parameters as in Ashes sim
    plt.plot(ash.span, ash.p_t[i], '-', label=f'Ashes {ash.windspeeds[i]}') # Plot the comparisons between Ashes and our BEM
    plt.plot(ash.BEM_rs, ash.BEM_p_t, '-', label=f'BEM {ash.windspeeds[i]}')
    plt.legend()
    plt.title(f"Comparison, BEM and Ashes at {ash.windspeeds[i]} [m/s]")
    plt.xlabel('Span [m]')
    plt.ylabel('Distributed tangential load [N/m]')
    plt.grid()
    plt.savefig(os.getcwd() + f"\\BEM-comparison-outputs\\comparison-{ash.windspeeds[i]}mps.png")
    plt.cla()


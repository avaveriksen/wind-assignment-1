'''
This file opens the raw data from the Ashes simulations of loads on the blades.
In the Ashes sensor panel, 'Thrust force, distr.' and 'Torque force distr.' are the ones
I've extracted. The raw output has a lot of crap from all the other sensors.
I wrote remove_header() function to shave off the upper header rows of the .txt output files.
The Extracted outputs from Ashes is stored in the 2D arrays p_n and p_t, starting with 5 m/s, then 9 and so forth.
The corresponding radii on the blade are found in vector 'span'.

It would be fairly time consuming to do it for a large range of pitch angles,
so I did it only at angle = 0.
'''

import pandas as pd
import matplotlib.pyplot as plt
import os

windspeeds = [5, 9, 11, 20] # m/s
p_n = [0] * 4 # [[p_n 5m/s], [p_n 9 m/s], ...]
p_t = [0] * 4 # same

#Blade span
span = [0, 2.64302, 5.37977, 8.20274, 11.1031, 14.071, 17.0953, 20.1641, 23.2647, 26.3837, 29.5076, 32.6228, 35.7156, 38.773, 41.7824, 44.732, 47.6111, 50.4099, 53.1201, 55.7344, 58.247, 60.6534, 62.9501, 65.1352, 67.2076, 69.1675, 71.0159, 72.7545, 74.386, 75.9133, 77.3402, 78.6705, 79.9085, 81.0585, 82.1252, 83.113, 84.0265, 84.8703, 85.6487, 86.366]

def remove_header(input_file, output_file):
    n_lines = 23;

    with open(input_file, "r") as f_in:
        lines = f_in.readlines()

    # Skip the first n_lines
    remaining = lines[n_lines:]

    with open(output_file, "w") as f_out:
        f_out.writelines(remaining)

files = ["5ms0deg.txt",
         "9ms0deg.txt",
         "11ms0deg.txt",
         "20ms0deg.txt"]
for file in files:
    file_name = os.getcwd() + "\\ashes-raw\\" + file
    out = os.getcwd() + "\\ashes-headers-removed\\" + file
    #remove_header(file_name, out)

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
    p_n[i] = df.iloc[n_rows - 1, end-(2*40):end - 40 + 1].values.tolist()
    p_t[i] = df.iloc[n_rows - 1, end-(40):end + 1].values.tolist()

#plt.plot(span, p_t[3])
#plt.show()


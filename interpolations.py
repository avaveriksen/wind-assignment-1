
import pandas as pd
import numpy as np

files = ["FFA-W3-241.txt", "FFA-W3-301.txt", "FFA-W3-360.txt", "FFA-W3-480.txt", "FFA-W3-600.txt", "cylinder.txt"]


tc1 = 0.301
tc2 = 0.360
tc_target = 0.350


tcs = [100, 86.05, ]

df1 = pd.read_csv(file1, delimiter='\t', header=None, names = ['alpha', 'cl', 'cd', 'cm'])
df2 = pd.read_csv(file2, delimiter='\t', header=None, names = ['alpha', 'cl', 'cd', 'cm'])

target = (tc1 + tc2) / 2;
df_interp = df1 + (target - tc1) * ((df2 - df1) / (tc2 - tc1))
df_interp["alpha"] = df1["alpha"]

df1.head
df2.head
df_interp.head




a = 0;
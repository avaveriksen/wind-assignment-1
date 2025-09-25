from windtools import Tools, Interpolations # Our toolbox of helper functions
import numpy as np 
import pandas as pd

V0=8.0 # m/s
R=89.166 # m
omega=0.6729 # degrees
fpath1="interpolated-tables/FFA_W3-1.0.csv"
beta=14.5 # degrees
c=5.38
theta_p=0 # degrees
r=2.8 # m
pn, pt, a, aprime, F = Tools.bem_single_element(r, c, beta, V0, omega, theta_p , fpath1, R)

print("pn is"+str(pn/1000))
print("pt is"+str(pt/1000))


for row in range(1, len(pd.read_csv("bladedat.txt", sep=r"\s+", header=None)) - 1):
                r, c, beta, tc, R, fpath = Tools.read_blade_data(row)
                pn, pt, a, aprime, F = Tools.bem_single_element(r, c, beta, V0, omega, theta_p, fpath, R)
                print("pn is"+str(pn/1000)+"     "+"pt is"+str(pt/1000))
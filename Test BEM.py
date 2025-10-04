from windtools import Tools # Our toolbox of helper functions
import numpy as np 
import pandas as pd
import scipy.integrate as integrate
from scipy.interpolate import RegularGridInterpolator

tc_values = [0.241, 0.301, 0.360, 0.480, 0.600, 1.0]  # Thickness-to-chord ratios
file_paths = [
        "FFA-W3-241.txt",
        "FFA-W3-301.txt",
        "FFA-W3-360.txt",
        "FFA-W3-480.txt",
        "FFA-W3-600.txt",
        "cylinder.txt",
        ]

    # Load them interpolation matrixs
cl_interp_func, cd_interp_func = Tools.load_airfoil_data(tc_values, file_paths)
print("✅ Airfoil data loaded and interpolation functions ready.")


V0=8.0 # m/s
R=89.166 # m
omega=0.6729 # rad/s
B=3 # number of blades
tc=0.241
beta=-3.3575 # degrees
c=1.1836
theta_p=0 # degree
r=88.4487 # m
Lambda = (omega*R)/V0

pn, pt, a, a_prime, F, alpha_deg, Cl, Cd= Tools.BEM(R,B,Lambda,theta_p,beta,c,r,V0,tc,cl_interp_func, cd_interp_func)
print("pn =", pn/1000, "N")
print("pt =", pt/1000, "N")










for row in range(1, len(pd.read_csv("bladedat.txt", sep=r"\s+", header=None)) - 1):
    r, c, beta, tc, R, fpath = Tools.read_blade_data(row)

    pn, pt, a, a_prime, F, alpha_deg, Cl, Cd = Tools.BEM(R, B, Lambda, theta_p, beta, c, r, V0, tc, cl_interp_func, cd_interp_func)
    print(f"r = {round(r, 3)} | pn = {round(pn / 1000, 3)} kN | pt = {round(pt / 1000, 3)} kN alpha"+str(alpha_deg)+"Cl="+str(Cl))

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import os

df = pd.read_csv("bladedat.txt", sep=r"\s+",header=None)

# Select a row (e.g. row 1 = second row in table)
Row = 1

# Assign variables
r = df.iloc[Row, 0]       # radius [m]
c = df.iloc[Row, 1]       # chord length [m]
beta = df.iloc[Row, 2]    # twist angle [deg]
tc = df.iloc[Row, 3] / 100  # convert % to fraction
R = df.iloc[-1,0] # Rotor radius [m]
# Example: fixed values you define, not from table
tip_speed_ratio = 6
theta_p = 0  # pitch angle [deg], if you want to set it yourself

V0=10 # This is an arbritrary value for the wind speed which we can normalize for later, when comparing optimal cp 
omega = tip_speed_ratio * V0 / R 

#Interpolation function for Cl and Cd based on alpha and t/c
def interpolate_from_table(tc_target, alpha, folder="interpolated-tables"): 
    # Format t/c to match filename (e.g., 0.241 -> 'FFA_W3-0.241.csv')
    rounded_tc = round(tc_target, 4) # Round to avoid floating point issues
    fname = f"FFA_W3-{rounded_tc}.csv"
    fpath = os.path.join(folder, fname)
    if not os.path.exists(fpath):
        raise FileNotFoundError(f"No table for t/c={tc_target} found in {folder}")
    df = pd.read_csv(fpath)
    # Interpolate for alpha
    cl_interp = interp1d(df['alpha'], df['cl'], kind='linear', fill_value='extrapolate')
    cd_interp = interp1d(df['alpha'], df['cd'], kind='linear', fill_value='extrapolate')
    cl = cl_interp(alpha)
    cd = cd_interp(alpha)
    return cl, cd

def bem_single_element(r, c, beta, tip_speed_ratio, theta_p, tc_target, R, B=3):
    rho = 1.225
    a = 0.0
    aprime = 0.0
    f_relax = 0.1
    tol = 1e-6
    max_iter = 100
    
    for i in range(max_iter):
        # Flow angle
        phi = np.arctan(((1 - a) ) / ((1 + aprime) * tip_speed_ratio * r/R))
        
        theta = np.deg2rad(beta + theta_p)  # Total pitch angle in radians
        # Angle of attack
        alpha = np.rad2deg(theta-phi)
        # Tip loss factor
        F = (2 / np.pi) * np.arccos(np.exp(-(B * (R - r)) / (2 * r * np.sin(abs(phi)))))
        F = max(F, 1e-5)
        
        #Cl and Cd from interpolation of t/c 
        Cl, Cd = interpolate_from_table(tc_target, alpha)
        # Aero coefficients
        Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
        Ct = Cl * np.sin(phi) - Cd * np.cos(phi)
        sigma = (c * B) / (2 * np.pi * r)  # Solidity
        
        # Glauert correction
        if a < 1/3:
            a_ast = ((sigma * Cn) * (1 - a)) / (4 * F * np.sin(phi)**2)
        else:
            dCT = (((1 - a)**2) * Cn * sigma) / (np.sin(phi)**2)
            a_ast = 0.246 * (dCT / F) + 0.0586 * (dCT / F)**2 + 0.0883 * (dCT / F)**3
        
        a_new = f_relax * a_ast + a * (1 - f_relax)
        
        aprime_ast = (((sigma * Ct) * (1 + a_ast)) / (4 * F * np.sin(phi) * np.cos(phi)))
        aprime_new = f_relax * aprime_ast + aprime * (1 - f_relax)
        
        # Convergence check
        if abs(a_new - a) < tol and abs(aprime_new - aprime) < tol:
            a, aprime = a_new, aprime_new
            break
        
        a, aprime = a_new, aprime_new
    
    # Relative wind speed
    Vrel = V0 * np.sqrt((1 - a)**2 + (tip_speed_ratio * r/R * (1 + aprime))**2)
    pn = 0.5 * rho * (Vrel**2) * c * Cn
    pt = 0.5 * rho * (Vrel**2) * c * Ct
    return pn, pt, a, aprime, F


# Print to check
print("r =", r, "m")
print("c ="
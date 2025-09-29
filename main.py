from windtools import Tools, Interpolations # Our toolbox of helper functions
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
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


if __name__ == '__main__':
    '''main loop'''

    

    
    # Don't run this always, comment out when we have out tables
    #Interpolations.interp_tc() # Interpolate tables across tc

    #pn, pt, a, aprime, F = tools.bem_single_element(PROVIDE PARAMETERS)
    """
    r=88.45
    c=1.18
    beta=-3.36
    V0= 4.0 # m/s
    omega=0.226116 # rad/s
    theta_p=-4.0 # degrees
    file = 'interpolated-tables\FFA_W3-0.241.csv'
    """ 
    #Tip speed ratio defined variables and pitch angle
    V0=8.0 # m/s
    R=89.166 # m
    B=3 #number of blades
    pitch_grid = np.arange(-4, 3,0.5)         # -4..3
    tsr_grid   = np.arange(5, 10,1)         # 5..10
    P_kW = np.zeros((len(tsr_grid), len(pitch_grid)))
    
    for i, lam in enumerate(tsr_grid):
        for j, theta_p in enumerate(pitch_grid):
            omega = lam * V0 / R
            rs, pts = [], []
            for row in range(1, len(pd.read_csv("bladedat.txt", sep=r"\s+", header=None)) - 1):
                r, c, beta, tc, R, fpath = Tools.read_blade_data(row)
                pn, pt, a, a_prime, F, alpha_deg, Cl, Cd = Tools.BEM(R, B, lam, theta_p, beta, c, r, V0, tc, cl_interp_func, cd_interp_func)
                rs.append(float(r)); pts.append(float(pt))
            rs = np.array(rs); pts = np.array(pts)
            P_kW[i, j] = B*np.trapezoid(pts * rs * omega, rs)   # power in W

    print("Power matrix (kW):")
    #print(P_kW)
    #print("Max power (kW):", np.max(P_kW), "is at", P_kW[np.unravel_index(np.argmax(P_kW), P_kW.shape)])

    # Calculate power coefficient CP
    Cp_matrix = P_kW / (0.5 * 1.225 * np.pi * R**2 * V0**3)  # Power coefficient
    max_cp = np.max(Cp_matrix)
    max_cp_index = np.unravel_index(np.argmax(Cp_matrix), Cp_matrix.shape)
    print("Maximum Cp:", max_cp, "at TSR =", tsr_grid[max_cp_index[0]], "and Pitch =", pitch_grid[max_cp_index[1]])
    print(Cp_matrix)

    T_kN = np.zeros((len(tsr_grid), len(pitch_grid)))
   
    for i, lam in enumerate(tsr_grid):
        for j, theta_p in enumerate(pitch_grid):
            omega = lam * V0 / R
            rs, fts = [], []
            for row in range(1, len(pd.read_csv("bladedat.txt", sep=r"\s+", header=None)) - 1):
                r, c, beta, tc, R, fpath = Tools.read_blade_data(row)
                pn, pt, a, a_prime, F, alpha_deg, Cl, Cd = Tools.BEM(R, B, lam, theta_p, beta, c, r, V0, tc, cl_interp_func, cd_interp_func)
                rs.append(float(r)); fts.append(float(pn))
            rs = np.array(rs); fts = np.array(fts)
            T_kN[i, j] = np.trapezoid(fts, rs)   # thrust in N

    # Calculate thrust coefficient CT
    Ct_matrix = T_kN / (0.5 * 1.225 * np.pi * R**2 * V0**2)  # Thrust coefficient
    print(Cp_matrix)
    
    
"""    # Create grids for plotting
TSR, Pitch = np.meshgrid(pitch_grid, tsr_grid)

plt.figure(figsize=(8,6))
cp_contour = plt.contourf(TSR, Pitch, Cp_matrix, levels=20, cmap="viridis")
plt.colorbar(cp_contour, label="$C_p$")
plt.xlabel("Pitch angle θp [deg]")
plt.ylabel("Tip-speed ratio λ")
plt.title("Contour plot of $C_p(λ, θ_p)$")
plt.show()
"""
    # Create grids for plotting
TSR, Pitch = np.meshgrid(pitch_grid, tsr_grid)

plt.figure(figsize=(8,6))
ct_contour = plt.contourf(TSR, Pitch, Ct_matrix, levels=20, cmap="viridis")
plt.colorbar(ct_contour, label="$C_p$")
plt.xlabel("Pitch angle θp [deg]")
plt.ylabel("Tip-speed ratio λ")
plt.title("Contour plot of $C_T(λ, θ_p)$")
plt.show()
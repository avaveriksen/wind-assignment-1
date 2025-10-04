from windtools import Tools, Interpolations # Our toolbox of helper functions
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
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

    pitch_grid = np.arange(-1, 2.1, 0.1 )         # -4..3
    tsr_grid   = np.arange(8,9.6 , 0.1)         # 5..10
    P_kW = np.zeros((len(tsr_grid), len(pitch_grid)))
    print(P_kW)
    for i, lam in enumerate(tsr_grid):
        for j, theta_p in enumerate(pitch_grid):
            omega = lam * V0 / R
            rs, pts = [], []
            for row in range(1, len(pd.read_csv("bladedat.txt", sep=r"\s+", header=None)) - 1):
                r, c, beta, tc, R, fpath = Tools.read_blade_data(row)
                pn, pt, a, aprime, F = Tools.bem_single_element(r, c, beta, V0, omega, theta_p, fpath, R)
                rs.append(float(r)); pts.append(float(pt))
            rs = np.array(rs); pts = np.array(pts)
            P_kW[i, j] = np.trapezoid(pts * rs * omega, rs) / 1e3  # power in kW

    print("Power matrix (kW):")
    print(P_kW)
    print("Max power (kW):", np.max(P_kW), "is at", P_kW[np.unravel_index(np.argmax(P_kW), P_kW.shape)])

    # Calculate power coefficient Cp
    Cp_matrix = (P_kW * 1e3) / (0.5 * 1.225 * np.pi * R**2 * V0**3)  # Power coefficient
    print(Cp_matrix)
    
    
    # Create grids for plotting
    TSR, Pitch = np.meshgrid(pitch_grid, tsr_grid)

    plt.figure(figsize=(8,6))
    cp_contour = plt.contourf(TSR, Pitch, Cp_matrix, levels=20, cmap="viridis")
    plt.colorbar(cp_contour, label="$C_p$")
    plt.xlabel("Pitch angle θp [deg]")
    plt.ylabel("Tip-speed ratio λ")
    plt.title("Contour plot of $C_p(λ, θ_p)$")
    plt.show()
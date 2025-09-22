from windtools import Tools, Interpolations # Our toolbox of helper functions
import numpy as np 
import pandas as pd
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
    V0=4.0 # m/s
    omega=0.226116 # rad/s
    theta_p=-4.0 # degrees

    rs, pts = [], []

    for row in range (1,len(pd.read_csv("bladedat.txt", sep=r"\s+",header=None))-1): # Loop over all rows except the last one (rotor radius)
        r,c,beta,tc,R,fpath = Tools.read_blade_data(row) # Read blade data from file, provide row number
        pn, pt, a, aprime, F = Tools.bem_single_element(r,c,beta,V0,omega,theta_p,fpath,R) # Run BEM on a single element, provide parameters within the function
        rs.append(float(r)); pts.append(float(pt))
        print(f"Results at r={r} m: pn={pn} N/m, pt={pt} N/m, a={a}, aprime={aprime}, F={F}")

    rs = np.array(rs); pts = np.array(pts)
    print("radius=",rs,"thrust=", pts)
    print(f"Torque: {np.trapz(pts*rs,rs)} Nm")
    print(f"Power: {np.trapz(pts*rs*omega,rs)/10**3} kW")
import numpy as np
import os
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
class Tools:
    

    def load_airfoil_data(tc_values, file_paths):
        """
        Load airfoil data for multiple thickness-to-chord ratios (tc_values)
        and create 2D interpolators for Cl and Cd as functions of (tc, alpha).

        Parameters:
            tc_values   : list of thickness-to-chord ratios corresponding to files
            file_paths  : list of file paths for each tc

        Returns:
            cl_interp_func : RegularGridInterpolator for Cl(tc, alpha)
            cd_interp_func : RegularGridInterpolator for Cd(tc, alpha)
        """
        all_alpha = []  # list to collect alpha grids
        all_cl = []     # Cl data
        all_cd = []     # Cd data

        # Loop over files (each tc)
        for f in file_paths:
            # Read whitespace-separated file with no header
            df = pd.read_csv(f, sep=r"\s+", header=None, engine='python')

            # Columns: alpha[deg], Cl, Cd, Cm
            alpha = df.iloc[:, 0].values
            cl = df.iloc[:, 1].values
            cd = df.iloc[:, 2].values

            # Wrap alpha into [-180, 180]
            alpha = ((alpha + 180) % 360) - 180

            # Sort alpha and remove duplicates
            alpha_sorted, unique_idx = np.unique(alpha, return_index=True)
            cl_sorted = cl[unique_idx]
            cd_sorted = cd[unique_idx]

            all_alpha.append(alpha_sorted)
            all_cl.append(cl_sorted)
            all_cd.append(cd_sorted)

        # Make sure all alpha grids are identical
        # Take the first as reference
        alpha_grid = all_alpha[0]
        for a in all_alpha[1:]:
            if not np.allclose(a, alpha_grid):
                raise ValueError("Alpha grids are not identical across files!")

        # Build Cl and Cd arrays (shape: tc x alpha)
        cl_array = np.array(all_cl)  # shape: (n_tc, n_alpha)
        cd_array = np.array(all_cd)  # shape: (n_tc, n_alpha)

        # Create 2D interpolators
        cl_interp_func = RegularGridInterpolator(
            (tc_values, alpha_grid),
            cl_array,
            bounds_error=False,
            fill_value=None  # returns NaN if out-of-bounds
        )

        cd_interp_func = RegularGridInterpolator(
            (tc_values, alpha_grid),
            cd_array,
            bounds_error=False,
            fill_value=None
        )

        return cl_interp_func, cd_interp_func



    

    

    

    def BEM(R, B, Lambda, theta_p, beta, c, r, V0, tc, cl_interp_func, cd_interp_func, max_iter=1000):
        """
        Single blade element BEM calculation with safety clamps and alpha wrapping.
        
        Parameters:
            R               : rotor radius [m]
            B               : number of blades
            Lambda          : tip speed ratio
            theta_p         : pitch angle [deg]
            beta            : twist angle at this section [deg]
            c               : chord at this section [m]
            r               : radial position of element [m]
            V0              : free-stream wind speed [m/s]
            tc              : thickness-to-chord ratio
            cl_interp_func  : RegularGridInterpolator for Cl(tc, alpha)
            cd_interp_func  : RegularGridInterpolator for Cd(tc, alpha)
            max_iter        : maximum iterations for a, a_prime
            
        Returns:
            pn, pt         : normal and tangential forces [N]
            a, a_prime     : axial and tangential induction factors
            F              : Prandtl tip-loss factor
            alpha_deg      : angle of attack [deg]
            Cl, Cd         : lift and drag coefficients
        """

        # Initialize induction factors
        a = 0.0
        a_new = 0.0
        a_prime = 0.0
        a_prime_new = 0.0

        f = 0.1  # relaxation factor
        rho = 1.225
        eps = 1e-6  # convergence tolerance

        # Loop for iterative solution
        for count in range(max_iter):
            a = a_new
            a_prime = a_prime_new

            sigma = (c * B) / (2 * np.pi * r)  # local solidity

            # Relative flow angle
            phi = np.arctan((1 - a) * R / ((1 + a_prime) * Lambda * r))
            phi = np.clip(phi, 1e-6, np.pi / 2 - 1e-6)  # prevent division by zero

            # Angle of attack
            alpha = phi - np.deg2rad(theta_p) - np.deg2rad(beta)
            alpha_deg = np.rad2deg(alpha)

            # Wrap alpha into [-180, 180] for interpolation
            #alpha_deg = ((alpha_deg + 180) % 360) - 180

            # Interpolate Cl and Cd
            Cl = cl_interp_func([[tc, alpha_deg]])[0]
            Cd = cd_interp_func([[tc, alpha_deg]])[0]

            # Compute normal and tangential force coefficients
            Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
            Ct = Cl * np.sin(phi) - Cd * np.cos(phi)

            # Prandtl tip-loss factor
            F = (2 / np.pi) * np.arccos(
                np.clip(np.exp(-(B / 2) * (R - r) / (r * np.sin(np.abs(phi)))), 0.0, 1.0)
            )
            F = np.clip(F, 0.0001, 1.0)  # safety clamp to prevent division by zero

            # Madsen formula for axial induction factor
            CT_star = ((1 - a)**2 * Cn * sigma) / (F * np.sin(phi)**2)
            CT_star = np.clip(CT_star, -1e3, 1e3)  # safety clamp

            a_new = 0.246 * CT_star + 0.0586 * CT_star**2 + 0.0883 * CT_star**3
            a_new = f*np.clip(a_new, 0.0, 0.95)+(1-f)*a  # physical limit

            # Tangential induction
            a_prime_star = sigma * Ct * (1 + a_prime) / (4 * F * np.sin(phi) * np.cos(phi))
            a_prime_star = np.clip(a_prime_star, -1e3, 1e3)
            a_prime_new = f * a_prime_star + (1 - f) * a_prime
            a_prime_new = np.clip(a_prime_new, -1.0, 1.0)  # physical limit

            # Safety: stop if NaN occurs
            if np.isnan(a_new) or np.isnan(a_prime_new):
                a_new = a
                a_prime_new = a_prime
                break

            # Check convergence
            if np.abs(a_new - a) < eps and np.abs(a_prime_new - a_prime) < eps:
                break

        # Relative velocity at the section
        omega = Lambda * V0 / R
        Vrel = np.sqrt((V0 * (1 - a))**2 + (omega * r * (1 + a_prime))**2)

        # Force per unit span
        # If r is very close to the tip, set forces to zero (class recommendation)
        if np.isclose(r, R, atol=1e-3):
            pn = 0.0
            pt = 0.0
        else:
            pn = 0.5 * rho * Vrel**2 * c * Cn
            pt = 0.5 * rho * Vrel**2 * c * Ct

        return pn, pt, a, a_prime, F, alpha_deg, Cl, Cd







    
    def read_blade_data(Row):
        df = pd.read_csv("bladedat.txt", sep=r"\s+",header=None)

        # Assign variables from blade data table 
        r = df.iloc[Row, 0]       # radius [m]
        c = df.iloc[Row, 1]       # chord length [m]
        beta = df.iloc[Row, 2]    # twist angle [deg]
        tc = df.iloc[Row, 3] / 100  # convert % to fraction
        R = df.iloc[-1,0] # Rotor radius [m]
        rounded_tc = round(tc, 4) # Round to avoid floating point issues
        fname = f"FFA_W3-{rounded_tc}.csv"
        fpath = os.path.join("interpolated-tables", fname)
        if not os.path.exists(fpath):
            raise FileNotFoundError(f"No table for t/c={tc} found in interpolated-tables")
        
        return r,c,beta,tc,R,fpath
    
    def bem_single_element(r, c, beta, V0, omega, theta_p, file, R, B=3):
        '''Compute the BEM algorithm for a single element'''

         # Load aerodynamic data from the file
        df = pd.read_csv(file)

        rho = 1.225
        theta = np.deg2rad(theta_p + beta)
        a = 0.0
        aprime = 0.0
        f_relax = 0.1
        tol = 1e-8
        max_iter = 300
        
        for i in range(max_iter):
            # Flow angle
            phi = np.arctan(((1 - a) * V0) / ((1 + aprime) * omega * r))

            # Tip loss factor
            F = (2 / np.pi) * np.arccos(np.exp(-(B * (R - r)) / (2 * r * np.sin(abs(phi)))))
            F = max(F, 1e-7)

            alpha = (phi - theta)
            alpha_deg = np.rad2deg(alpha)

            # Interpolate Cl and Cd from the dataframe
            cl_interp = np.interp(alpha_deg, df['alpha'], df['cl'])
            cd_interp = np.interp(alpha_deg, df['alpha'], df['cd'])
            Cl = cl_interp
            Cd = cd_interp



            # Aero coefficients
            Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
            Ct = Cl * np.sin(phi) - Cd * np.cos(phi)
            sigma = (c * B) / (2 * np.pi * r)  # Solidity

            # Glauert correction
            if a < 1 / 3:
                a_ast = ((sigma * Cn) * (1 - a)) / (4 * F * np.sin(phi) ** 2)
            else:
                dCT = (((1 - a) ** 2) * Cn * sigma) / (np.sin(phi) ** 2)
                a_ast = 0.246 * (dCT / F) + 0.0586 * (dCT / F) ** 2 + 0.0883 * (dCT / F) ** 3

            a_new = f_relax * a_ast + a * (1 - f_relax)

            aprime_ast = (((sigma * Ct) * (1 + aprime)) / (4 * F * np.sin(phi) * np.cos(phi)))
            aprime_new = f_relax * aprime_ast + aprime * (1 - f_relax)

            # Convergence check
            if abs(a_new - a) < tol and abs(aprime_new - aprime) < tol:
                a, aprime = a_new, aprime_new
                break
            
            a, aprime = a_new, aprime_new
            

        # Relative wind speed
        Vrel = np.sqrt((V0 * (1 - a)) ** 2 + (omega * r * (1 + aprime)) ** 2)
        pn = 0.5 * rho * (Vrel ** 2) * c * Cn
        pt = 0.5 * rho * (Vrel ** 2) * c * Ct
        return pn, pt, a, aprime, F,alpha_deg, Cl   
        
    def extend_bladedata():
        source = pd.read_csv("bladedat.txt", sep=r"\s+", header=None)
        mylist = []
        for row in range(0, 18):
            fname = f"FFA_W3-{round(source.iloc[row, 3] / 100, 4)}.csv"
            mylist.append(fname)

        source[4] = mylist
        file_name = os.getcwd() + "\\extended_bladedata.csv"
        source.to_csv(file_name, index=False)
    
    
        

class Interpolations:
    def interp_tc(self):
        # data provided
        files = ["FFA-W3-241.txt",
                 "FFA-W3-301.txt",
                 "FFA-W3-360.txt",
                 "FFA-W3-480.txt",
                 "FFA-W3-600.txt",
                 "cylinder.txt"]
        # accompanying tc
        tc_list = [0.241,
                   0.301,
                   0.360,
                   0.480,
                   0.600,
                   1]
        # including tc's in 10MW blade
        tc_extended = [0.241,
                       0.2426,
                       0.2532,
                       0.2781,
                       0.3242,
                       0.4304,
                       0.611,
                       0.8605,
                       1]
        # data in one list [(tc, dataframe), (tc, dataframe)...]
        airfoil_data = []

        for i, tc in enumerate(tc_list):
            # fill airfoil_data with (tc, dataframe) tuples
            df = pd.read_csv(files[i], delimiter='\t', header=None, names=['alpha', 'cl', 'cd', 'cm'])
            airfoil_data.append((tc, df))

        for i, tc in enumerate(tc_extended):
            # The tc's for which were not provided tabels are linearly interpolated
            if not (tc in tc_list):
                # variables for linear slope calculation
                tc1 = airfoil_data[i - 1][0]
                df1 = airfoil_data[i - 1][1]
                tc2 = airfoil_data[i][0]
                df2 = airfoil_data[i][1]

                # linear interpolation
                df_interp = df1 + (tc - tc1) * ((df2 - df1) / (tc2 - tc1))

                # round values in dataframe
                df_interp = df_interp.__round__(4)

                # chop up, insert and concatenate airfoil_data with interp
                arr = airfoil_data[:i] + [(tc, df_interp)] + airfoil_data[i:]
                airfoil_data = arr

        for tup in airfoil_data:
            # save to csv files
            file_name = os.getcwd() + "\\interpolated-tables\\FFA_W3-" + str(round(tup[0], 4)) + ".csv"
            tup[1].to_csv(file_name, index=False)

import numpy as np
import os
import pandas as pd

class Tools:

    def bem_single_element(r, c, beta, V0, omega, theta_p, file, R=89.17, B=3):
        '''Compute the BEM algorithm for a single element'''

         # Load aerodynamic data from the file
        df = pd.read_csv(file)

        rho = 1.225
        theta = np.deg2rad(theta_p + beta)
        a = 0.0
        aprime = 0.0
        f_relax = 0.1
        tol = 1e-6
        max_iter = 100

        for i in range(max_iter):
            # Flow angle
            phi = np.arctan(((1 - a) * V0) / ((1 + aprime) * omega * r))

            # Tip loss factor
            F = (2 / np.pi) * np.arccos(np.exp(-(B * (R - r)) / (2 * r * np.sin(abs(phi)))))
            F = max(F, 1e-5)

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

            aprime_ast = (((sigma * Ct) * (1 + a_ast)) / (4 * F * np.sin(phi) * np.cos(phi)))
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
        return pn, pt, a, aprime, F

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

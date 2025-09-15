import numpy as np



def bem_single_element(r, c, beta, V0, omega, theta_p, Cl, Cd, R=89.17, B=3):
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
    Vrel = np.sqrt((V0 * (1 - a))**2 + (omega * r * (1 + aprime))**2)
    pn = 0.5 * rho * (Vrel**2) * c * Cn
    pt = 0.5 * rho * (Vrel**2) * c * Ct
    return pn, pt, a, aprime, F

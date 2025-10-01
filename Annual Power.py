from windtools import Tools, Interpolations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =============================================================
# --- Load airfoil data (kept from your original code) ---
# =============================================================
tc_values = [0.241, 0.301, 0.360, 0.480, 0.600, 1.0]
file_paths = [
    "FFA-W3-241.txt",
    "FFA-W3-301.txt",
    "FFA-W3-360.txt",
    "FFA-W3-480.txt",
    "FFA-W3-600.txt",
    "cylinder.txt",
]

cl_interp_func, cd_interp_func = Tools.load_airfoil_data(tc_values, file_paths)
print("✅ Airfoil data loaded and interpolation functions ready.")

# =============================================================
# --- Hidden main aerodynamic loop ---
# =============================================================

# Uncomment this block if you need to regenerate the power curve later
R = 89.166  # m
B = 3
pitch_grid = np.arange(-4, 3, 0.5)
tsr_grid = np.arange(5, 10, 1)
V_range = np.arange(3, 26, 1)
"""
V_list, Pmax_list = [], []
for V0 in V_range:
    print(f"Evaluating wind speed {V0:.1f} m/s...")
    P_kW = np.zeros((len(tsr_grid), len(pitch_grid)))

    for i, lam in enumerate(tsr_grid):
        for j, theta_p in enumerate(pitch_grid):
            omega = lam * V0 / R
            rs, pts = [], []
            for row in range(1, len(pd.read_csv('bladedat.txt', sep=r'\s+', header=None)) - 1):
                r, c, beta, tc, R, fpath = Tools.read_blade_data(row)
                pn, pt, a, a_prime, F, alpha_deg, Cl, Cd = Tools.BEM(
                    R, B, lam, theta_p, beta, c, r, V0, tc, cl_interp_func, cd_interp_func
                )
                rs.append(float(r))
                pts.append(float(pt))
            rs = np.array(rs)
            pts = np.array(pts)
            P_kW[i, j] = B * np.trapezoid(pts * rs * omega, rs) / 1000

    P_max_kW = np.max(P_kW)
    V_list.append(V0)
    Pmax_list.append(P_max_kW)
    print(f"→ Max power at {V0:.1f} m/s: {P_max_kW:.2f} kW")
P_array = np.array(Pmax_list) / 1000  # Convert to MW

# =============================================================
# --- Stored test arrays (from your power curve summary) ---
# =============================================================
"""
V_array = np.arange(3, 26, 1)  # 3–25 m/s
P_array = np.array([
    195.49, 463.38, 905.03, 1563.89, 2483.41, 3707.01, 5278.14, 7240.25,
    9636.77, 12511.15, 15906.83, 19867.25, 24435.84, 29656.06, 35571.35, 42225.14,
    49660.88, 57922.00, 67051.96, 77094.18, 88092.12, 100089.22, 113128.91
]) / 1000  # Convert to MW

# =============================================================
# --- Weibull distribution + AEP calculation ---
# =============================================================

# --- 2. Weibull PDF function
def weibull_pdf(V, A, k):
    """Compute Weibull probability density function."""
    return (k / A) * (V / A)**(k - 1) * np.exp(-(V / A)**k)

# --- 3. Site-specific Weibull parameters
A = 9.0   # scale parameter [m/s]
k = 1.9   # shape parameter [-]

# --- 4. Compute probability density for each wind speed
f_V = weibull_pdf(V_array, A, k)

# --- 5. Rated power and capped curve
P_rated_MW = 10.64
P_array_capped = np.minimum(P_array, P_rated_MW)

# --- 6. Compute expected (mean) power output by integrating P(V)*f(V)
P_expected_uncapped = np.trapz(P_array * f_V, V_array)
P_expected_capped   = np.trapz(P_array_capped * f_V, V_array)

# --- 7. Compute AEPs
hours_per_year = 8760
AEP_uncapped_MWh = P_expected_uncapped * hours_per_year
AEP_capped_MWh   = P_expected_capped * hours_per_year

# --- 8. Compute probability and energy beyond 20 m/s
mask_20 = V_array <= 20
AEP20_capped_MWh = np.trapz(P_array_capped[mask_20] * f_V[mask_20], V_array[mask_20]) * hours_per_year
prob_above_20 = np.trapz(f_V[V_array > 20], V_array[V_array > 20]) * 100  # in %

# --- 9. Print results
print("\n=== Energy Production Summary ===")
print(f"Expected mean power (uncapped): {P_expected_uncapped:.3f} MW")
print(f"Expected mean power (capped @ {P_rated_MW:.2f} MW): {P_expected_capped:.3f} MW")
print(f"AEP (uncapped): {AEP_uncapped_MWh:.2f} MWh/year")
print(f"AEP (capped):   {AEP_capped_MWh:.2f} MWh/year")
print(f"AEP (capped, ≤20 m/s): {AEP20_capped_MWh:.2f} MWh/year")
print(f"Probability of V > 20 m/s: {prob_above_20:.3f}%")
print(f"Energy loss by capping (vs uncapped): {(AEP_uncapped_MWh - AEP_capped_MWh):.2f} MWh ({100*(AEP_uncapped_MWh - AEP_capped_MWh)/AEP_uncapped_MWh:.2f}%)")

# --- 10. Plot Weibull + both power curves
fig, ax1 = plt.subplots(figsize=(9,6))

# Bars for Weibull PDF
bars = ax1.bar(V_array, f_V, width=0.8, color='skyblue', alpha=0.8, label='Weibull PDF')
for i, v in enumerate(V_array):
    if v > 20:
        bars[i].set_color('orange')

ax1.set_xlabel("Wind speed V [m/s]")
ax1.set_ylabel("Probability density f(V)")
ax1.set_ylim(0, max(f_V)*1.3)
ax1.legend(loc="upper left")

# Power curves on secondary axis
ax2 = ax1.twinx()
ax2.plot(V_array, P_array, 'o-', color='navy', label="Power curve (uncapped)")
ax2.plot(V_array, P_array_capped, '--', color='red', linewidth=2, label=f"Capped @ {P_rated_MW} MW")
ax2.set_ylabel("Power (MW)")
ax2.legend(loc="upper right")

plt.title("Weibull Distribution and Power Curve\n(orange bars = V > 20 m/s)")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
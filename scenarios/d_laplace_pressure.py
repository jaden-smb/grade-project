#!/usr/bin/env python3

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import time
import numpy as np
import matplotlib.pyplot as plt

from lbm.core import LBMSimulator

OUTPUT_DIR = "output/scenario_d_laplace"


def _bulk_density(rho_arr, rho_min, rho_max, phase='liquid', bulk_fraction=0.10):

    rng = rho_max - rho_min
    if phase == 'liquid':
        mask = rho_arr > rho_max - bulk_fraction * rng
    else:
        mask = rho_arr < rho_min + bulk_fraction * rng
    if not mask.any():
        # Fall back to single-cell extreme if the bulk region is empty
        return float(rho_arr.max()) if phase == 'liquid' else float(rho_arr.min())
    return float(rho_arr[mask].mean())


def run(tau=1.0, G=-5.0, rho_liquid=2.0, rho_gas=0.1,
        radii=None, nx=200, ny=200, num_steps=15000):

    if radii is None:
        radii = [20, 30, 40, 50, 60]

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}/")
    print(f"Grid: {nx}x{ny}, tau={tau}, G={G}, steps={num_steps}")
    print(f"Radii: {radii}")

    cs2 = 1.0 / 3.0

    def psi(rho):
        return 1.0 - np.exp(-1.5 * rho)

    def sc_pressure(rho_bulk):
        """Full SC EoS: p = rho*cs² + G*cs²/2 * psi(rho)²"""
        return rho_bulk * cs2 + G * cs2 / 2.0 * psi(rho_bulk) ** 2

    inv_R_list  = []
    delta_p_list = []
    rows = []  # for the summary table

    for R in radii:
        print(f"\n  Radius R={R} ...")
        lbm = LBMSimulator(nx, ny, tau, G, rho_liquid, rho_gas)
        lbm.initialize_droplet(nx // 2, ny // 2, float(R))

        t_start = time.perf_counter()
        for _ in range(num_steps):
            lbm.step()
        elapsed = time.perf_counter() - t_start

        rho_arr = np.array(lbm.get_density()).reshape((ny, nx))
        rho_max = float(rho_arr.max())
        rho_min = float(rho_arr.min())

        # --- Effective radius from liquid area ---
        rho_threshold = (rho_max + rho_min) / 2.0
        area = float((rho_arr > rho_threshold).sum())
        R_eff = float(np.sqrt(area / np.pi)) if area > 0 else float(R)

        # --- Bulk-averaged densities (excludes the interface) ---
        rho_liq_bulk = _bulk_density(rho_arr, rho_min, rho_max, phase='liquid')
        rho_gas_bulk = _bulk_density(rho_arr, rho_min, rho_max, phase='gas')

        # --- Pressures from the full SC EoS ---
        p_liq   = sc_pressure(rho_liq_bulk)
        p_gas   = sc_pressure(rho_gas_bulk)
        delta_p = p_liq - p_gas

        print(f"    rho_liq_bulk={rho_liq_bulk:.4f}  rho_gas_bulk={rho_gas_bulk:.4f}")
        print(f"    R_eff={R_eff:.1f} lu  Δp={delta_p:.5f}  t={elapsed:.1f}s")

        if R_eff > 0:
            inv_R_list.append(1.0 / R_eff)
            delta_p_list.append(delta_p)
            rows.append((R, rho_liq_bulk, rho_gas_bulk, R_eff, delta_p,
                         1.0 / R_eff, elapsed))

    if len(inv_R_list) < 2:
        print("Not enough data points for Laplace fit.")
        return

    inv_R   = np.array(inv_R_list)
    delta_p = np.array(delta_p_list)

    # --- Linear fit: Δp = σ * (1/R) + intercept ---
    coeffs    = np.polyfit(inv_R, delta_p, 1)
    sigma     = coeffs[0]
    intercept = coeffs[1]

    # --- R² ---
    delta_p_pred = np.polyval(coeffs, inv_R)
    ss_res    = float(np.sum((delta_p - delta_p_pred) ** 2))
    ss_tot    = float(np.sum((delta_p - delta_p.mean()) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')

    print(f"\n  Fitted surface tension  σ         = {sigma:.4f} lu·[pressure]")
    print(f"  Intercept (should be ~0)           = {intercept:.5f}")
    print(f"  R²                                 = {r_squared:.5f}")

    # --- Summary table ---
    print("\n  R0(lu)  rho_liq   rho_gas   Reff(lu)  Δp        1/Reff    t(s)")
    for row in rows:
        print(f"  {row[0]:6d}  {row[1]:.4f}    {row[2]:.4f}    "
              f"{row[3]:6.1f}    {row[4]:+.5f}  {row[5]:.5f}   {row[6]:.1f}")

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(inv_R, delta_p, color='b', s=70, zorder=3, label='Simulated Δp')
    x_fit = np.linspace(0, inv_R.max() * 1.15, 200)
    ax.plot(x_fit, np.polyval(coeffs, x_fit), 'r--', linewidth=1.5,
            label=f'Linear fit:  σ = {sigma:.4f},  R² = {r_squared:.4f}')
    ax.set_xlabel('1 / R$_\\mathrm{eff}$  (lu$^{-1}$)', fontsize=12)
    ax.set_ylabel('$\\Delta p = p_\\mathrm{liq} - p_\\mathrm{gas}$', fontsize=12)
    ax.set_title('Laplace Pressure Test\n(Young–Laplace: $\\Delta p = \\sigma / R$)',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # Annotate intercept as a note inside the plot
    ax.annotate(f'Intercept = {intercept:.4f}',
                xy=(0.05, 0.12), xycoords='axes fraction', fontsize=9,
                color='gray')

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/laplace_test.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved: {OUTPUT_DIR}/laplace_test.png")

    return inv_R_list, delta_p_list, sigma, r_squared


if __name__ == "__main__":
    run()

#!/usr/bin/env python3
"""Scenario D — Laplace Pressure Test.

Simulates droplets of five different radii and fits Δp vs 1/R to extract the
surface tension σ, validating the Young-Laplace relation (Δp = σ/R).
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import time
import numpy as np
import matplotlib.pyplot as plt

from lbm.core import LBMSimulator

OUTPUT_DIR = "output/scenario_d_laplace"


def run(tau=1.0, G=-5.0, rho_liquid=2.0, rho_gas=0.1,
        radii=None, nx=200, ny=200, num_steps=5000):

    if radii is None:
        radii = [20, 30, 40, 50, 60]

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}/")
    print(f"Grid: {nx}x{ny}, tau={tau}, G={G}, steps={num_steps}")
    print(f"Radii: {radii}")

    inv_R_list = []
    delta_p_list = []

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

        rho_threshold = (rho_max + rho_min) / 2.0
        area = float((rho_arr > rho_threshold).sum())
        R_eff = float(np.sqrt(area / np.pi)) if area > 0 else float(R)

        cs2 = 1.0 / 3.0
        def psi(rho):
            return 1.0 - np.exp(-1.5 * rho)
        p_liq = rho_max * cs2 + G * cs2 / 2.0 * psi(rho_max)**2
        p_gas = rho_min * cs2 + G * cs2 / 2.0 * psi(rho_min)**2
        delta_p = p_liq - p_gas

        print(f"    rho_liq={rho_max:.4f}, rho_gas={rho_min:.4f}, "
              f"R_eff={R_eff:.1f}, Δp={delta_p:.5f}, t={elapsed:.1f}s")

        if R_eff > 0:
            inv_R_list.append(1.0 / R_eff)
            delta_p_list.append(delta_p)

    if len(inv_R_list) < 2:
        print("Not enough data points for Laplace fit.")
        return

    inv_R = np.array(inv_R_list)
    delta_p = np.array(delta_p_list)

    coeffs = np.polyfit(inv_R, delta_p, 1)
    sigma = coeffs[0]
    print(f"\n  Fitted surface tension σ = {sigma:.4f} (lattice units)")

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(inv_R, delta_p, color='b', s=60, zorder=3, label='Simulated Δp')
    x_fit = np.linspace(0, inv_R.max() * 1.1, 100)
    ax.plot(x_fit, np.polyval(coeffs, x_fit), 'r--', linewidth=1.5,
            label=f'Linear fit  σ = {sigma:.4f}')
    ax.set_xlabel('1 / R_eff  (lu⁻¹)', fontsize=12)
    ax.set_ylabel('Δp = p_liq − p_gas', fontsize=12)
    ax.set_title('Laplace Pressure Test\n(Young–Laplace: Δp = σ/R)',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/laplace_test.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {OUTPUT_DIR}/laplace_test.png")

    return inv_R_list, delta_p_list, sigma


if __name__ == "__main__":
    run()

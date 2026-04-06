#!/usr/bin/env python3

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import time
import numpy as np
import matplotlib.pyplot as plt

from lbm.core import LBMSimulator

OUTPUT_DIR = "output/scenario_c_coexistence"


def run(nx=150, ny=150, tau=1.0, radius=35, num_steps=3000,
        rho_liquid=2.0, rho_gas=0.1,
        G_values=None):

    if G_values is None:
        G_values = [-4.0, -4.5, -5.0, -5.5, -6.0]

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}/")
    print(f"Grid: {nx}x{ny}, tau={tau}, radius={radius}, steps={num_steps}")
    print(f"G values: {G_values}")

    eq_rho_liquid = []
    eq_rho_gas = []

    for G in G_values:
        print(f"\n  Running G = {G} ...")
        lbm = LBMSimulator(nx, ny, tau, G, rho_liquid, rho_gas)
        lbm.initialize_droplet(nx // 2, ny // 2, radius)

        t_start = time.perf_counter()
        for _ in range(num_steps):
            lbm.step()
        elapsed = time.perf_counter() - t_start

        rho_final = np.array(lbm.get_density()).reshape((ny, nx))
        eq_rho_liquid.append(float(rho_final.max()))
        eq_rho_gas.append(float(rho_final.min()))
        print(f"    rho_liq={rho_final.max():.4f}, rho_gas={rho_final.min():.4f}, "
              f"t={elapsed:.2f}s ({nx*ny*num_steps/elapsed/1e6:.2f} MLUPS)")

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(G_values, eq_rho_liquid, 'b-o', linewidth=1.5, markersize=7, label='Liquid phase (max ρ)')
    ax.plot(G_values, eq_rho_gas, 'r-o', linewidth=1.5, markersize=7, label='Gas phase (min ρ)')
    ax.set_xlabel('G (cohesion parameter)', fontsize=12)
    ax.set_ylabel('Equilibrium Density', fontsize=12)
    ax.set_title('Coexistence Curve: Liquid and Gas Densities vs. G',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/coexistence_curve.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved: {OUTPUT_DIR}/coexistence_curve.png")

    return G_values, eq_rho_liquid, eq_rho_gas


if __name__ == "__main__":
    run()

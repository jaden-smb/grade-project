#!/usr/bin/env python3
"""Scenario B — G-Ramp Evaporation Analogy.

Ramps the cohesion parameter G from -5.0 to -3.7 over 8000 steps.
Past the critical point (Gc ≈ -4.0) the droplet dissolves, simulating evaporation.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from lbm.runner import run_simulation
from lbm.visualization import plot_density_profile, plot_velocity_field

OUTPUT_DIR = "output/scenario_b_evaporation"


def run():
    lbm, rho, elapsed = run_simulation(
        nx=200, ny=200,
        tau=1.0,
        G=-5.0, G_final=-3.7,
        rho_liquid=2.0, rho_gas=0.1,
        radius=40,
        num_steps=8000,
        save_plots=True,
        animate=True,
        view_surface=True,
        output_dir=OUTPUT_DIR,
    )
    plot_density_profile(rho, 200, 200, OUTPUT_DIR)
    plot_velocity_field(lbm, 200, 200, OUTPUT_DIR)


if __name__ == "__main__":
    run()

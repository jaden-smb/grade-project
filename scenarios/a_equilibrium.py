#!/usr/bin/env python3

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from lbm.runner import run_simulation
from lbm.visualization import plot_density_profile, plot_velocity_field

OUTPUT_DIR = "output/scenario_a_equilibrium"


def run():
    lbm, rho, elapsed = run_simulation(
        nx=200, ny=200,
        tau=1.0,
        G=-5.0,
        rho_liquid=2.0, rho_gas=0.1,
        radius=40,
        num_steps=5000,
        save_plots=True,
        animate=True,
        view_surface=True,
        export_obj=True,
        output_dir=OUTPUT_DIR,
    )
    plot_density_profile(rho, 200, 200, OUTPUT_DIR)
    plot_velocity_field(lbm, 200, 200, OUTPUT_DIR)


if __name__ == "__main__":
    run()

#!/usr/bin/env python3
"""Entry point for the LBM Shan-Chen two-phase simulation.

Usage:
  python main.py                  # run all four scenarios
  python main.py --scenario a     # steady-state equilibrium
  python main.py --scenario b     # G-ramp evaporation analogy
  python main.py --scenario c     # coexistence curve sweep
  python main.py --scenario d     # Laplace pressure test
"""

import argparse
import importlib
import sys
import os

# Make src/ importable for lbm package, and root importable for scenarios package
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

SCENARIOS = {
    "a": ("scenarios.a_equilibrium",     "Scenario A — Steady-State Droplet Equilibrium"),
    "b": ("scenarios.b_evaporation",     "Scenario B — G-Ramp Evaporation Analogy"),
    "c": ("scenarios.c_coexistence_curve", "Scenario C — Parameter Sweep for Coexistence Curve"),
    "d": ("scenarios.d_laplace_pressure", "Scenario D — Laplace Pressure Test"),
}


def main():
    parser = argparse.ArgumentParser(
        description="LBM Shan-Chen Two-Phase Simulation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Scenarios:\n" + "\n".join(f"  {k}  {v[1]}" for k, v in SCENARIOS.items()),
    )
    parser.add_argument(
        "--scenario", "-s",
        choices=[*SCENARIOS.keys(), "all"],
        default="all",
        metavar="{a,b,c,d,all}",
        help="Scenario to run (default: all)",
    )
    args = parser.parse_args()

    to_run = list(SCENARIOS.keys()) if args.scenario == "all" else [args.scenario]

    for key in to_run:
        module_path, label = SCENARIOS[key]
        print(f"\n{'=' * 60}")
        print(label)
        print("=" * 60)
        mod = importlib.import_module(module_path)
        mod.run()

    print("\nAll scenarios complete.")


if __name__ == "__main__":
    main()

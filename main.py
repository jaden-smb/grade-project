#!/usr/bin/env python3

import argparse
import importlib
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

SCENARIOS = {
    "a": ("scenarios.a_equilibrium",       "Scenario A — Steady-State Droplet Equilibrium"),
    "b": ("scenarios.b_evaporation",       "Scenario B — G-Ramp Evaporation Analogy"),
    "c": ("scenarios.c_coexistence_curve", "Scenario C — Coexistence Curve + Maxwell Construction"),
    "d": ("scenarios.d_laplace_pressure",  "Scenario D — Laplace Pressure Test"),
    "e": ("scenarios.e_eos_comparison",    "Scenario E — EoS Comparison (SC vs. Carnahan-Starling)"),
    "f": ("scenarios.f_parameter_study",   "Scenario F — Systematic Parameter Study"),
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
        metavar="{a,b,c,d,e,f,all}",
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

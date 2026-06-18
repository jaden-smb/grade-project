#!/usr/bin/env python3

import argparse
import importlib
import sys
import os
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))


class _Tee:
    def __init__(self, stream, path):
        self._stream = stream
        self._file = open(path, "w", encoding="utf-8")

    def write(self, data):
        self._stream.write(data)
        self._file.write(data)

    def flush(self):
        self._stream.flush()
        self._file.flush()

    def close(self):
        self._file.close()

    def __getattr__(self, name):
        return getattr(self._stream, name)

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

    log_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = os.path.join(log_dir, f"simulation_{timestamp}.txt")

    tee = _Tee(sys.stdout, log_path)
    sys.stdout = tee
    try:
        print(f"Log: {log_path}")
        print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

        for key in to_run:
            module_path, label = SCENARIOS[key]
            print(f"\n{'=' * 60}")
            print(label)
            print("=" * 60)
            mod = importlib.import_module(module_path)
            mod.run()

        print("\nAll scenarios complete.")
        print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    finally:
        sys.stdout = tee._stream
        tee.close()
        print(f"Log saved to: {log_path}")


if __name__ == "__main__":
    main()

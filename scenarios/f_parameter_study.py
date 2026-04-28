#!/usr/bin/env python3
"""Scenario F — Systematic Parameter Study.

Documents simulation behaviour under different parameter choices
(thesis objective: "Documentación del comportamiento bajo diferentes parámetros").

Two sweeps:
  1. tau sweep at fixed G=-5.0, varying tau in {0.7, 0.8, 1.0, 1.2, 1.5}
     Measures: equilibrium density contrast, effective droplet radius,
               max spurious velocity, MLUPS throughput.

  2. Initial density contrast sweep at fixed G=-5.0, tau=1.0:
     rho_liquid in {1.5, 2.0, 2.5, 3.0} at fixed rho_gas=0.1.

Results are printed as a table and saved as a multi-panel figure.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from lbm.core import LBMSimulator

OUTPUT_DIR = "output/scenario_f_parameter_study"

# Common simulation settings
_NX, _NY = 120, 120
_RADIUS = 28
_NUM_STEPS = 2000
_G_BASE = -5.0
_RHO_GAS_BASE = 0.1
_RHO_LIQ_BASE = 2.0


def _effective_radius(rho_2d, rho_threshold):
    area = float((rho_2d > rho_threshold).sum())
    return float(np.sqrt(area / np.pi)) if area > 0 else 0.0


def _run_one(nx, ny, tau, G, rho_liquid, rho_gas, radius, num_steps):
    """Return dict of measured metrics after running one simulation."""
    lbm = LBMSimulator(nx, ny, tau, G, rho_liquid, rho_gas)
    lbm.initialize_droplet(nx // 2, ny // 2, radius)

    t0 = time.perf_counter()
    for _ in range(num_steps):
        lbm.step()
    elapsed = time.perf_counter() - t0

    rho = np.array(lbm.get_density()).reshape((ny, nx))
    ux  = np.array(lbm.get_velocity_x()).reshape((ny, nx))
    uy  = np.array(lbm.get_velocity_y()).reshape((ny, nx))

    rho_max = float(rho.max())
    rho_min = float(rho.min())
    contrast = rho_max / rho_min if rho_min > 1e-6 else float('nan')

    rho_threshold = (rho_max + rho_min) / 2.0
    r_eff = _effective_radius(rho, rho_threshold)

    speed = np.sqrt(ux ** 2 + uy ** 2)
    max_spur_vel = float(speed.max())

    # Mean spurious velocity in the interface region (20–80% of density range)
    rho_range = rho_max - rho_min
    if rho_range > 1e-6:
        interface_mask = (rho > rho_min + 0.2 * rho_range) & \
                         (rho < rho_min + 0.8 * rho_range)
        mean_spur_vel = float(speed[interface_mask].mean()) if interface_mask.any() else 0.0
    else:
        mean_spur_vel = 0.0

    mlups = nx * ny * num_steps / elapsed / 1e6

    # Detect instability: negative density means distribution functions blew up;
    # max spurious velocity > 1.0 (Mach ~ 1.7) also signals breakdown.
    stable = rho_min >= 0.0 and max_spur_vel <= 1.0

    return {
        'rho_max':        rho_max,
        'rho_min':        rho_min,
        'contrast':       contrast,
        'r_eff':          r_eff,
        'max_spur_vel':   max_spur_vel,
        'mean_spur_vel':  mean_spur_vel,
        'mlups':          mlups,
        'stable':         stable,
    }


def _print_table(title, header_row, rows, col_fmts):
    """Print a simple text table.

    String values (e.g. 'UNSTABLE') are printed as-is without applying the
    format spec, so callers can mark bad rows without special-casing the helper.
    """
    col_w = [max(len(h), 10) for h in header_row]
    sep = "  ".join("-" * w for w in col_w)
    hdr = "  ".join(f"{h:>{w}}" for h, w in zip(header_row, col_w))
    print(f"\n{title}")
    print("=" * len(sep))
    print(hdr)
    print(sep)
    for row in rows:
        fmt_row = [
            v if isinstance(v, str) else fmt.format(v)
            for v, fmt in zip(row, col_fmts)
        ]
        print("  ".join(f"{s:>{w}}" for s, w in zip(fmt_row, col_w)))
    print("=" * len(sep))


def run(nx=_NX, ny=_NY, num_steps=_NUM_STEPS):

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}/")
    print(f"Grid: {nx}x{ny}, steps per run: {num_steps}")

    # ── Sweep 1: tau ────────────────────────────────────────────────────────
    tau_values = [0.7, 0.8, 1.0, 1.2, 1.5]
    print(f"\n{'─'*60}")
    print(f"Sweep 1: tau ∈ {tau_values}  (G={_G_BASE}, rho_liq={_RHO_LIQ_BASE})")
    print(f"{'─'*60}")

    sweep1 = []
    for tau in tau_values:
        print(f"  tau={tau:.1f} ...", end=" ", flush=True)
        m = _run_one(nx, ny, tau, _G_BASE, _RHO_LIQ_BASE,
                     _RHO_GAS_BASE, _RADIUS, num_steps)
        sweep1.append(m)
        if m['stable']:
            print(f"contrast={m['contrast']:.2f}  r_eff={m['r_eff']:.1f}  "
                  f"u_max={m['max_spur_vel']:.2e}  {m['mlups']:.1f} MLUPS")
        else:
            print(f"UNSTABLE  rho_min={m['rho_min']:.4f}  "
                  f"u_max={m['max_spur_vel']:.2e}  {m['mlups']:.1f} MLUPS")

    _UNSTABLE = 'UNSTABLE'
    _print_table(
        "Sweep 1 — tau variation (G=-5.0, rho_liq=2.0, rho_gas=0.1)",
        ["tau", "rho_liq", "rho_gas", "contrast", "r_eff", "u_max(spur)", "MLUPS"],
        [
            [tau, m['rho_max'], m['rho_min'], m['contrast'],
             m['r_eff'], m['max_spur_vel'], m['mlups']]
            if m['stable'] else
            [tau, _UNSTABLE, _UNSTABLE, _UNSTABLE,
             _UNSTABLE, _UNSTABLE, m['mlups']]
            for tau, m in zip(tau_values, sweep1)
        ],
        ["{:.1f}", "{:.4f}", "{:.4f}", "{:.2f}", "{:.1f}", "{:.2e}", "{:.1f}"],
    )

    unstable_taus = [tau for tau, m in zip(tau_values, sweep1) if not m['stable']]
    if unstable_taus:
        tau_str = ' and '.join(f'tau={t}' for t in unstable_taus)
        print(f"\nNote: {tau_str} are numerically unstable at G={_G_BASE}"
              f" (negative densities or max |u| > 1.0).")
        min_stable = min(t for t, m in zip(tau_values, sweep1) if m['stable'])
        print(f"      Stable range: tau >= ~{min_stable} at this G value.")

    # ── Sweep 2: initial rho_liquid ─────────────────────────────────────────
    rho_liq_values = [1.5, 2.0, 2.5, 3.0]
    print(f"\n{'─'*60}")
    print(f"Sweep 2: rho_liquid ∈ {rho_liq_values}  (G={_G_BASE}, tau=1.0)")
    print(f"{'─'*60}")

    sweep2 = []
    for rho_liq in rho_liq_values:
        print(f"  rho_liquid={rho_liq:.1f} ...", end=" ", flush=True)
        m = _run_one(nx, ny, 1.0, _G_BASE, rho_liq,
                     _RHO_GAS_BASE, _RADIUS, num_steps)
        sweep2.append(m)
        print(f"contrast={m['contrast']:.2f}  r_eff={m['r_eff']:.1f}  "
              f"u_max={m['max_spur_vel']:.2e}  {m['mlups']:.1f} MLUPS")

    _print_table(
        "Sweep 2 — initial rho_liquid variation (G=-5.0, tau=1.0, rho_gas=0.1)",
        ["rho_liq0", "rho_liq_eq", "rho_gas_eq", "contrast", "r_eff", "u_max(spur)", "MLUPS"],
        [[rho_liq, m['rho_max'], m['rho_min'], m['contrast'],
          m['r_eff'], m['max_spur_vel'], m['mlups']]
         for rho_liq, m in zip(rho_liq_values, sweep2)],
        ["{:.1f}", "{:.4f}", "{:.4f}", "{:.2f}", "{:.1f}", "{:.2e}", "{:.1f}"],
    )

    # ── Figure ───────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 4, figsize=(18, 9))

    def _plot_sweep(row, x_vals, metrics, x_label, sweep_title):
        keys   = ['contrast', 'r_eff', 'max_spur_vel', 'mlups']
        titles = ['Density Contrast ρ_liq/ρ_gas',
                  'Effective Radius (lu)',
                  'Max Spurious Velocity (lu/step)',
                  'Throughput (MLUPS)']
        colors = ['steelblue', 'darkorange', 'crimson', 'seagreen']

        stable_mask = [m.get('stable', True) for m in metrics]
        x_stable    = [x for x, s in zip(x_vals, stable_mask) if s]
        x_unstable  = [x for x, s in zip(x_vals, stable_mask) if not s]

        for col, (key, title, color) in enumerate(zip(keys, titles, colors)):
            ax = axes[row][col]
            vals_stable = [m[key] for m, s in zip(metrics, stable_mask) if s]
            ax.plot(x_stable, vals_stable, color=color,
                    marker='o', linewidth=1.5, markersize=7)
            # Mark unstable x positions with vertical dotted lines
            for ux in x_unstable:
                ax.axvline(x=ux, color='red', linestyle=':', alpha=0.6, linewidth=1.2)
            ax.set_xlabel(x_label, fontsize=10)
            ax.set_title(title, fontsize=10, fontweight='bold')
            ax.grid(True, alpha=0.3)
            if col == 0:
                ax.set_ylabel(sweep_title, fontsize=9, labelpad=2)
            # Annotation on first panel only
            if col == 0 and x_unstable:
                unstable_str = ', '.join(f'{x}' for x in x_unstable)
                ax.text(0.03, 0.97,
                        f'Not plotted (unstable):\n{x_label} = {unstable_str}',
                        transform=ax.transAxes, fontsize=7, va='top', color='red',
                        bbox=dict(boxstyle='round,pad=0.3',
                                  facecolor='mistyrose', alpha=0.8))

    _plot_sweep(0, tau_values, sweep1, 'τ', 'Sweep 1')
    _plot_sweep(1, rho_liq_values, sweep2, 'ρ_liquid', 'Sweep 2')

    plt.suptitle(
        'Parameter Study — Effect of τ and Initial Density Contrast\n'
        f'(G = {_G_BASE},  {nx}×{ny} grid,  {num_steps} steps)',
        fontsize=13, fontweight='bold',
    )
    plt.tight_layout()
    outfile = f'{OUTPUT_DIR}/parameter_study.png'
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved: {outfile}")

    return {'sweep_tau': sweep1, 'sweep_rho_liq': sweep2,
            'tau_values': tau_values, 'rho_liq_values': rho_liq_values}


if __name__ == "__main__":
    run()

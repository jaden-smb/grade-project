#!/usr/bin/env python3
"""Scenario E — EoS Comparison (Objective 1: Shan-Chen variants).

Runs the same droplet simulation with two pseudopotential formulations:
  EoS 0: Original Shan-Chen   psi = 1 - exp(-1.5*rho)
  EoS 1: Carnahan-Starling    psi = sqrt(2*(p_CS - rho*cs2) / (G*cs2))
               (Yuan & Schaefer 2006)

For each EoS the coexistence curve (equilibrium rho_liq and rho_gas vs G)
is measured and plotted on the same axes together with the density contrast
ratio.  An inset table summarises the contrast at each G value.
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

OUTPUT_DIR = "output/scenario_e_eos_comparison"

# Carnahan-Starling critical temperature (a=0.5, b=2.0, R=1):  Tc ≈ 0.09433
_CS_TC = 0.09433
_CS_T_DEFAULT = 0.7 * _CS_TC   # well into two-phase region


def _run_one(nx, ny, tau, G, rho_liquid, rho_gas, radius, num_steps,
             eos_type, T_cs):
    """Run a single simulation and return (rho_liq, rho_gas, MLUPS)."""
    lbm = LBMSimulator(nx, ny, tau, G, rho_liquid, rho_gas,
                       eos_type=eos_type, T=T_cs)
    lbm.initialize_droplet(nx // 2, ny // 2, radius)

    t0 = time.perf_counter()
    for _ in range(num_steps):
        lbm.step()
    elapsed = time.perf_counter() - t0

    rho = np.array(lbm.get_density()).reshape((ny, nx))
    mlups = nx * ny * num_steps / elapsed / 1e6
    return float(rho.max()), float(rho.min()), mlups


def run(nx=120, ny=120, tau=1.0, radius=28, num_steps=3000,
        rho_liquid=2.0, rho_gas=0.1,
        G_values=None, T_cs=_CS_T_DEFAULT):

    if G_values is None:
        G_values = [-4.5, -5.0, -5.5, -6.0]

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}/")
    print(f"Grid: {nx}x{ny}, tau={tau}, radius={radius}, steps={num_steps}")
    print(f"G values: {G_values}")
    print(f"CS temperature T = {T_cs:.5f}  (Tc = {_CS_TC:.5f},  T/Tc = {T_cs/_CS_TC:.3f})")

    EOS_LABELS = {0: "Original Shan-Chen", 1: "Carnahan-Starling (Yuan-Schaefer)"}
    results = {}  # eos_type -> {'G': [], 'rho_liq': [], 'rho_gas': [], 'contrast': []}

    for eos_type in [0, 1]:
        label = EOS_LABELS[eos_type]
        print(f"\n{'─'*55}")
        print(f"  EoS {eos_type}: {label}")
        print(f"{'─'*55}")
        rho_liqs, rho_gases, contrasts = [], [], []

        for G in G_values:
            print(f"  G = {G:+.1f} ...", end=" ", flush=True)
            rl, rg, mlups = _run_one(nx, ny, tau, G, rho_liquid, rho_gas,
                                     radius, num_steps, eos_type, T_cs)
            contrast = rl / rg if rg > 1e-6 else float('nan')
            rho_liqs.append(rl)
            rho_gases.append(rg)
            contrasts.append(contrast)
            print(f"rho_liq={rl:.4f}  rho_gas={rg:.4f}  "
                  f"contrast={contrast:.2f}  ({mlups:.1f} MLUPS)")

        results[eos_type] = {
            'G': G_values,
            'rho_liq': rho_liqs,
            'rho_gas': rho_gases,
            'contrast': contrasts,
        }

    # ── Plot ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    colors = {0: ('#1f77b4', '#d62728'),   # blue/red for SC
              1: ('#2ca02c', '#ff7f0e')}    # green/orange for CS
    markers = {0: 'o', 1: 's'}

    # Left panel: coexistence curve
    ax = axes[0]
    for eos_type, label in EOS_LABELS.items():
        r = results[eos_type]
        cliq, cgas = colors[eos_type]
        mk = markers[eos_type]
        ax.plot(r['G'], r['rho_liq'], color=cliq, marker=mk, linewidth=1.5,
                markersize=7, label=f'{label} — liquid')
        ax.plot(r['G'], r['rho_gas'], color=cgas, marker=mk, linewidth=1.5,
                linestyle='--', markersize=7, label=f'{label} — gas')
    ax.set_xlabel('G (cohesion parameter)', fontsize=12)
    ax.set_ylabel('Equilibrium Density', fontsize=12)
    ax.set_title('Coexistence Curve by EoS', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Right panel: density contrast ratio
    ax = axes[1]
    for eos_type, label in EOS_LABELS.items():
        r = results[eos_type]
        cliq, _ = colors[eos_type]
        mk = markers[eos_type]
        ax.plot(r['G'], r['contrast'], color=cliq, marker=mk, linewidth=1.5,
                markersize=7, label=label)
    ax.set_xlabel('G (cohesion parameter)', fontsize=12)
    ax.set_ylabel('Density contrast ratio  ρ_liq / ρ_gas', fontsize=12)
    ax.set_title('Density Contrast by EoS', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.suptitle(
        'EoS Comparison: Original Shan-Chen vs. Carnahan-Starling\n'
        f'(CS: T/Tc = {T_cs/_CS_TC:.2f},  a = 0.5,  b = 2.0)',
        fontsize=13, fontweight='bold', y=1.01,
    )
    plt.tight_layout()
    outfile = f'{OUTPUT_DIR}/eos_comparison.png'
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved: {outfile}")

    # Print summary table
    print("\n" + "=" * 65)
    print(f"{'Density Contrast Summary':^65}")
    print("=" * 65)
    header = f"{'G':>6}  {'SC contrast':>12}  {'CS contrast':>12}  {'Δ contrast':>12}"
    print(header)
    print("-" * 65)
    sc = results[0]
    cs = results[1]
    for i, G in enumerate(G_values):
        sc_c = sc['contrast'][i]
        cs_c = cs['contrast'][i]
        delta = sc_c - cs_c
        print(f"{G:>6.1f}  {sc_c:>12.2f}  {cs_c:>12.2f}  {delta:>+12.2f}")
    print("=" * 65)
    print("Note: With the Yuan-Schaefer formulation, CS coexistence densities")
    print("depend only on temperature (T/Tc), not on G — hence the flat CS curves.")
    print("The original SC EoS density contrast varies with G because psi encodes")
    print("the EoS implicitly through the exponential form.")

    return results


if __name__ == "__main__":
    run()

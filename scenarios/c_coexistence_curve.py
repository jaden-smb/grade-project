#!/usr/bin/env python3

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import quad

from lbm.core import LBMSimulator

OUTPUT_DIR = "output/scenario_c_coexistence"


def _p_sc(rho, G):
    """Shan-Chen equation of state: p = rho*cs2 + G*cs2/2 * psi(rho)^2."""
    cs2 = 1.0 / 3.0
    psi = 1.0 - np.exp(-1.5 * rho)
    return rho * cs2 + G * cs2 / 2.0 * psi ** 2


def maxwell_construction(G, rho_max=4.0, n_grid=8000):
    """Return (rho_gas, rho_liq) from the Maxwell equal-area rule for the SC EoS.

    Solves simultaneously:
      (1)  p(rho_g) = p(rho_l)                     [equal pressure]
      (2)  ∫ from rho_g to rho_l (p(ρ) - p_coex) / ρ² dρ = 0
                                                    [equal chemical potential]

    Equation (2) is the equal-area rule in the (1/ρ, p) plane — the correct
    thermodynamic Maxwell construction for a barotropic EoS.  Using
    ∫(p − p_coex) dρ = 0 (equal area in the (ρ, p) plane) is WRONG and
    converges to the trivial degenerate root ρ_g = ρ_l.

    Returns (None, None) if no coexistence is found (|G| too small).
    """
    rho_grid = np.linspace(1e-3, rho_max, n_grid)
    p_grid = _p_sc(rho_grid, G)

    # Locate dp/dρ < 0 (spinodal) region
    dp = np.diff(p_grid) / np.diff(rho_grid)
    unstable = np.where(dp < 0)[0]
    if len(unstable) == 0:
        return None, None

    rho_sp_lo = rho_grid[unstable[0]]
    rho_sp_hi = rho_grid[unstable[-1] + 1]

    # Use initial guesses well inside the stable monotonic branches —
    # far from each other to avoid the trivial root ρ_g = ρ_l.
    rg0 = min(0.02, rho_sp_lo * 0.15)          # deep in gas branch
    rl0 = max(2.5, rho_sp_hi * 1.8)            # deep in liquid branch
    rl0 = min(rl0, rho_max * 0.92)             # stay within grid

    def equations(vars):
        rg, rl = vars
        if rg <= 0 or rl <= rg:
            return [1e6, 1e6]
        pg = float(_p_sc(rg, G))
        pl = float(_p_sc(rl, G))
        eq_pressure = pg - pl
        p_coex = pg
        # Maxwell equal-area in the (1/ρ, p) plane: ∫(p − p_coex)/ρ² dρ = 0
        integral, _ = quad(
            lambda r: (float(_p_sc(r, G)) - p_coex) / (r * r),
            rg, rl, limit=200,
        )
        return [eq_pressure, integral]

    try:
        sol, info, ier, _ = fsolve(equations, [rg0, rl0], full_output=True)
        rg, rl = sol
        residual = np.max(np.abs(info['fvec']))
        # Reject if fsolve did not converge, residuals too large,
        # densities unphysical, or trivial root (rho_liq / rho_gas < 2).
        if (ier != 1 or residual > 1e-8 or rg <= 0 or rl <= rg
                or rl > rho_max or rl / rg < 2.0):
            return None, None
        return float(rg), float(rl)
    except Exception:
        return None, None


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

    # Compute Maxwell construction analytical predictions
    print("\n  Computing Maxwell construction (analytical) ...")
    mx_rho_gas = []
    mx_rho_liq = []
    for G in G_values:
        rg, rl = maxwell_construction(G)
        mx_rho_gas.append(rg)
        mx_rho_liq.append(rl)
        if rg is not None:
            # Verification: print residuals so convergence is transparent.
            p_coex_v = _p_sc(rg, G)
            dp_check = _p_sc(rl, G) - p_coex_v
            area_check, _ = quad(
                lambda r: (_p_sc(r, G) - p_coex_v) / (r * r), rg, rl,
                limit=200,
            )
            print(f"    G={G:+.1f}: rho_gas={rg:.4f}, rho_liq={rl:.4f}"
                  f"  [verify: Δp={dp_check:.2e}, area={area_check:.2e}]")
        else:
            print(f"    G={G:+.1f}: no coexistence found")

    fig, ax = plt.subplots(figsize=(9, 6))

    # Simulated data (solid lines + filled markers)
    ax.plot(G_values, eq_rho_liquid, 'b-o', linewidth=1.5, markersize=7,
            label='Simulated — liquid phase')
    ax.plot(G_values, eq_rho_gas, 'r-o', linewidth=1.5, markersize=7,
            label='Simulated — gas phase')

    # Maxwell construction (dashed lines + open markers)
    mx_G_valid = [G for G, rg in zip(G_values, mx_rho_gas) if rg is not None]
    mx_rl_valid = [rl for rl in mx_rho_liq if rl is not None]
    mx_rg_valid = [rg for rg in mx_rho_gas if rg is not None]

    if mx_G_valid:
        ax.plot(mx_G_valid, mx_rl_valid, 'b--s', linewidth=1.2, markersize=7,
                markerfacecolor='none', label='Maxwell construction — liquid')
        ax.plot(mx_G_valid, mx_rg_valid, 'r--s', linewidth=1.2, markersize=7,
                markerfacecolor='none', label='Maxwell construction — gas')

    ax.set_xlabel('G (cohesion parameter)', fontsize=12)
    ax.set_ylabel('Equilibrium Density', fontsize=12)
    ax.set_title('Coexistence Curve: Simulated vs. Maxwell Construction (original SC EoS)',
                 fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/coexistence_curve.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved: {OUTPUT_DIR}/coexistence_curve.png")

    return G_values, eq_rho_liquid, eq_rho_gas


if __name__ == "__main__":
    run()

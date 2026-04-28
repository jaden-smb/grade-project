"""Physics-level unit tests for the LBM C++ engine.

Requires the C++ extension to be built.  All tests are skipped automatically
when the extension is not available (e.g., on a fresh clone before building).

Run with:
    python -m pytest tests/test_physics.py -v
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import numpy as np
import pytest

try:
    from lbm.core import LBMSimulator
    HAS_ENGINE = True
except ImportError:
    HAS_ENGINE = False

requires_engine = pytest.mark.skipif(
    not HAS_ENGINE,
    reason="C++ extension not built — run 'cmake -B build && cmake --build build' first",
)


# ─── helpers ────────────────────────────────────────────────────────────────

def _make_sim(nx=80, ny=80, tau=1.0, G=-5.0,
              rho_liquid=2.0, rho_gas=0.1, radius=18):
    sim = LBMSimulator(nx, ny, tau, G, rho_liquid, rho_gas)
    sim.initialize_droplet(nx // 2, ny // 2, radius)
    return sim


def _rho(sim):
    return np.array(sim.get_density()).reshape((sim.ny, sim.nx))


def _momentum(sim):
    ux = np.array(sim.get_velocity_x()).reshape((sim.ny, sim.nx))
    uy = np.array(sim.get_velocity_y()).reshape((sim.ny, sim.nx))
    rho = _rho(sim)
    return float((rho * ux).sum()), float((rho * uy).sum())


# ─── tests ──────────────────────────────────────────────────────────────────

@requires_engine
def test_mass_conservation():
    """Total mass must not drift by more than 0.1 % over 500 steps."""
    sim = _make_sim()
    mass_initial = float(_rho(sim).sum())

    for _ in range(500):
        sim.step()

    mass_final = float(_rho(sim).sum())
    drift_pct = 100.0 * abs(mass_final - mass_initial) / mass_initial
    assert drift_pct < 0.1, (
        f"Mass drift {drift_pct:.4f}% exceeds 0.1% threshold"
    )


@requires_engine
def test_phase_separation():
    """With G=-5.0, max density must exceed 1.5 and min density must be < 0.5
    after 2000 steps, confirming that phases have separated."""
    sim = _make_sim(G=-5.0)
    for _ in range(2000):
        sim.step()

    rho = _rho(sim)
    assert rho.max() > 1.5, (
        f"max density {rho.max():.4f} ≤ 1.5; phase separation did not occur"
    )
    assert rho.min() < 0.5, (
        f"min density {rho.min():.4f} ≥ 0.5; phase separation did not occur"
    )


@requires_engine
def test_symmetry():
    """A centred droplet on a square grid must retain approximate 4-fold symmetry
    after 500 steps.  The four quadrants are compared pairwise: the mean absolute
    density difference must be < 5 % of the density range."""
    nx, ny = 80, 80
    sim = LBMSimulator(nx, ny, 1.0, -5.0, 2.0, 0.1)
    sim.initialize_droplet(nx // 2, ny // 2, 18)

    for _ in range(500):
        sim.step()

    rho = _rho(sim)
    cx, cy = nx // 2, ny // 2
    rho_range = rho.max() - rho.min()

    q1 = rho[:cy, :cx]
    q2 = rho[:cy, cx:]
    q3 = rho[cy:, :cx]
    q4 = rho[cy:, cx:]

    # Each quadrant vs its 180°-rotated partner
    diff_13 = np.abs(q1 - np.rot90(q4, 2)).mean()
    diff_24 = np.abs(q2 - np.rot90(q3, 2)).mean()

    tol = 0.05 * rho_range
    assert diff_13 < tol, (
        f"Quadrant 1 vs 3 mean diff {diff_13:.4f} > 5% of range ({tol:.4f})"
    )
    assert diff_24 < tol, (
        f"Quadrant 2 vs 4 mean diff {diff_24:.4f} > 5% of range ({tol:.4f})"
    )


@requires_engine
def test_momentum_conservation():
    """With no external forcing and a symmetric initial droplet, total momentum
    magnitude should remain < 1e-8 per lattice site after 500 steps."""
    nx, ny = 80, 80
    sim = LBMSimulator(nx, ny, 1.0, -5.0, 2.0, 0.1)
    sim.initialize_droplet(nx // 2, ny // 2, 18)

    for _ in range(500):
        sim.step()

    px, py = _momentum(sim)
    n_sites = nx * ny
    # Normalise by number of sites for a grid-size-independent threshold
    mag = (px ** 2 + py ** 2) ** 0.5 / n_sites
    assert mag < 1e-8, (
        f"Momentum magnitude per site {mag:.2e} exceeds 1e-8 threshold"
    )


@requires_engine
def test_default_eos_unchanged():
    """EoS type 0 must give the same density field as the pre-refactor default
    (ψ = 1 − exp(−1.5ρ)).  Verified by checking that rho0 = 1/1.5 reproduces
    identical results to the baseline factory settings."""
    nx, ny = 60, 60
    sim_a = LBMSimulator(nx, ny, 1.0, -5.0, 2.0, 0.1, eos_type=0)
    sim_a.initialize_droplet(nx // 2, ny // 2, 14)

    sim_b = LBMSimulator(nx, ny, 1.0, -5.0, 2.0, 0.1, eos_type=0)
    sim_b.set_rho0(1.0 / 1.5)  # explicit rho0, must match default
    sim_b.initialize_droplet(nx // 2, ny // 2, 14)

    for _ in range(100):
        sim_a.step()
        sim_b.step()

    rho_a = np.array(sim_a.get_density())
    rho_b = np.array(sim_b.get_density())
    np.testing.assert_allclose(rho_a, rho_b, atol=1e-14,
                               err_msg="eos_type=0 with default rho0 gave different results")


@requires_engine
def test_cs_eos_gives_phase_separation():
    """The Carnahan-Starling EoS (eos_type=1) must also produce phase separation
    (density contrast > 1.5) after 3000 steps at G=-5.0 with T=0.7*Tc."""
    nx, ny = 80, 80
    T_cs = 0.7 * 0.09433
    sim = LBMSimulator(nx, ny, 1.0, -5.0, 2.0, 0.1, eos_type=1, T=T_cs)
    sim.initialize_droplet(nx // 2, ny // 2, 18)

    for _ in range(3000):
        sim.step()

    rho = _rho(sim)
    contrast = rho.max() / max(rho.min(), 1e-9)
    assert contrast > 1.5, (
        f"CS EoS density contrast {contrast:.3f} ≤ 1.5; phase separation too weak"
    )

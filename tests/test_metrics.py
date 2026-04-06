"""Unit tests for lbm.metrics — no C++ extension required."""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from lbm.metrics import compute_metrics


def test_perfect_circle_metrics():
    """A circular liquid region should have circularity close to 1 and correct radius."""
    ny, nx = 100, 100
    rho = np.zeros((ny, nx))
    cy, cx, R = 50, 50, 20

    ys, xs = np.ogrid[:ny, :nx]
    mask = (xs - cx)**2 + (ys - cy)**2 <= R**2
    rho[mask] = 2.0
    rho[~mask] = 0.1

    threshold = 1.05
    total_mass, radius, max_density, aspect_ratio, circularity = compute_metrics(rho, threshold)

    expected_radius = np.sqrt(mask.sum() / np.pi)
    assert abs(radius - expected_radius) < 1.0, f"radius {radius:.2f} != expected {expected_radius:.2f}"
    assert circularity > 0.85, f"circularity {circularity:.3f} should be > 0.85 for a circle"
    assert abs(max_density - 2.0) < 1e-9
    assert abs(total_mass - (2.0 * mask.sum() + 0.1 * (~mask).sum())) < 1e-6


def test_empty_field_returns_zero_radius():
    """When no pixels exceed the threshold, radius and circularity should be 0."""
    rho = np.full((50, 50), 0.1)
    total_mass, radius, max_density, aspect_ratio, circularity = compute_metrics(rho, threshold=1.0)

    assert radius == 0.0
    assert circularity == 0.0
    assert aspect_ratio == 0.0


def test_mass_conservation():
    """Total mass is always the sum of all density values regardless of threshold."""
    rng = np.random.default_rng(42)
    rho = rng.uniform(0.05, 2.5, size=(80, 80))
    total_mass, *_ = compute_metrics(rho, threshold=1.0)
    assert abs(total_mass - float(rho.sum())) < 1e-6


def test_aspect_ratio_circular():
    """A perfect circle should have aspect ratio close to 1."""
    ny, nx = 100, 100
    rho = np.zeros((ny, nx))
    ys, xs = np.ogrid[:ny, :nx]
    mask = (xs - 50)**2 + (ys - 50)**2 <= 25**2
    rho[mask] = 2.0

    _, _, _, aspect_ratio, _ = compute_metrics(rho, threshold=1.0)
    assert aspect_ratio > 0.9, f"aspect ratio {aspect_ratio:.3f} should be close to 1 for a circle"


if __name__ == "__main__":
    test_perfect_circle_metrics()
    test_empty_field_returns_zero_radius()
    test_mass_conservation()
    test_aspect_ratio_circular()
    print("All tests passed.")

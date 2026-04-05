# Quick Start Guide

## Prerequisites Check

Before building, ensure you have:
```bash
# Check Python version (need 3.6+)
python3 --version

# Check CMake version (need 3.12+)
cmake --version

# Create and activate virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install Python dependencies in venv
pip install numpy matplotlib pybind11
```

## Build and Run (CMake Method)

```bash
# Activate venv first (if using one)
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Option 1: Use the build script (automatically activates venv if present)
./build.sh

# Option 2: Manual build
mkdir -p build && cd build
cmake ..
cmake --build .
# Module is automatically placed in project root by CMake config
cd ..

# Run simulation (make sure venv is activated)
python python/simulate.py
```

## Build and Run (setup.py Method)

```bash
# Build extension
python setup.py build_ext --inplace

# Run simulation
python python/simulate.py
```

## Expected Output

The simulation runs four scenarios sequentially:

1. **Scenario A** (`output/scenario_a_equilibrium/`): 200×200 grid, G=-5.0, 5000 steps. Steady-state droplet equilibrium — demonstrates spontaneous phase separation into a stable liquid droplet in gas.
2. **Scenario B** (`output/scenario_b_evaporation/`): 200×200 grid, G ramped -5.0→-3.7, 8000 steps. Evaporation analogy — droplet dissolves as cohesion is ramped past the critical point (G_c ≈ -4.0).
3. **Scenario C** (`output/scenario_c_coexistence/`): 150×150 grid, 3000 steps per G value. Sweeps G ∈ {-4.0, -4.5, -5.0, -5.5, -6.0} and plots the coexistence curve.
4. **Scenario D** (`output/scenario_d_laplace/`): 200×200 grid, 5000 steps per radius. Tests five droplet radii (20–60 lu) and fits Δp vs 1/R using the full Shan-Chen EoS to extract surface tension σ.

Each scenario saves density plots, surface view images, animated GIFs, and a `metrics.png` with radius, mass, max density, aspect ratio, and circularity over time. Performance (MLUPS) is printed to the console at the end of each run.

## Quick Parameter Tuning

Edit `python/simulate.py`, function `main()`:

- **Faster simulation**: Reduce `nx`, `ny` (e.g., 100x100)
- **More detail**: Increase `nx`, `ny` (e.g., 400x400)
- **Stronger phase separation**: Make `G` more negative (e.g., -6.5)
- **Weaker phase separation**: Make `G` less negative (e.g., -4.5)
- **Higher viscosity**: Reduce `tau` (e.g., 0.7)
- **Lower viscosity**: Increase `tau` (e.g., 1.3)

## Troubleshooting

**Import Error**: Make sure the `.so` (Linux/Mac) or `.pyd` (Windows) file is in the project root.

**Compilation Error**: Check that you have C++17 support and PyBind11 installed.

**Unstable Simulation**: Try reducing `tau` or making `G` less negative.


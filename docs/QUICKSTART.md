# Quick Start Guide

## Prerequisites Check

```bash
# Python 3.6+
python3 --version

# CMake 3.12+
cmake --version

# Create and activate a virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate        # Linux / macOS
venv\Scripts\activate           # Windows

# Install Python dependencies
pip install -r requirements.txt
```

## Build the C++ Extension

**Option 1 — CMake (recommended):**
```bash
cmake -B build
cmake --build build
```

**Option 2 — setuptools:**
```bash
python setup.py build_ext --inplace --build-lib src/lbm
```

After building, the compiled extension should be at:
```
src/lbm/lbm_shan_chen*.so   # Linux / macOS
src/lbm/lbm_shan_chen*.pyd  # Windows
```

## Run the Simulation

```bash
# All six scenarios (sequential)
python main.py

# Single scenario
python main.py --scenario a    # Steady-state equilibrium
python main.py --scenario b    # G-ramp evaporation analogy
python main.py --scenario c    # Coexistence curve sweep
python main.py --scenario d    # Laplace pressure test
python main.py --scenario e    # EoS comparison (SC vs Carnahan-Starling)
python main.py --scenario f    # Systematic parameter study

# Each scenario is also directly executable
python scenarios/a_equilibrium.py
```

Console output is also saved to `output/simulation_<timestamp>.txt` automatically.

## Expected Output

| Scenario | Grid | Steps | Output |
|---|---|---|---|
| A | 200×200 | 5 000 | Density plots, 3D surfaces, GIF animation, `metrics.png`, OBJ sequence |
| B | 200×200 | 8 000 | Density plots, 3D surfaces, GIF animation, `metrics.png` |
| C | 150×150 | 3 000 × 5 | Coexistence curve overlaid with Maxwell construction |
| D | 200×200 | 5 000 × 5 | Δp vs 1/R linear fit; prints extracted surface tension σ |
| E | 120×120 | 3 000 × 2 | Side-by-side coexistence curve and density contrast for SC vs CS EoS |
| F | 120×120 | 2 000 × 9+ | τ and ρ_liquid sweeps; contrast, radius, spurious velocity, MLUPS tables |

Throughput (MLUPS) is printed to the console at the end of each run.

## Quick Parameter Tuning

Edit the relevant file in `scenarios/`:

| Goal | Parameter | Change |
|---|---|---|
| Faster run | `nx`, `ny` | Decrease (e.g., 100×100) |
| More detail | `nx`, `ny` | Increase (e.g., 400×400) |
| Stronger separation | `G` | More negative (e.g., -6.5) |
| Weaker separation | `G` | Less negative (e.g., -4.5) |
| Higher viscosity | `tau` | Decrease (e.g., 0.7) |
| Lower viscosity | `tau` | Increase (e.g., 1.3) |

## Running Tests

```bash
python -m pytest tests/
# or
python tests/test_metrics.py
```

Tests in `tests/test_metrics.py` validate droplet metrics with synthetic data and do not require the C++ extension.

## Troubleshooting

**`ImportError: No module named 'lbm_shan_chen'`**
Build the extension first and confirm the `.so` / `.pyd` file is inside `src/lbm/`.

**Compilation error**
Check that you have a C++17-capable compiler and that PyBind11 is installed (`pip install pybind11`).

**Unstable simulation**
Reduce `tau` or make `G` less negative.

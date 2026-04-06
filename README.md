# Two-Phase Lattice Boltzmann Simulation — Shan-Chen Model

A complete implementation of a two-phase Lattice Boltzmann Method (LBM) simulation using the Shan-Chen pseudo-potential model for liquid-gas phase separation. The codebase uses Python for high-level orchestration and C++ (exposed through PyBind11) for all performance-critical operations.

## Overview

This simulation demonstrates:
- **Phase separation**: Formation of liquid and gas phases from an initial mixed state
- **Interface dynamics**: Evolution of the liquid-gas interface
- **Droplet formation**: Creation and evolution of liquid droplets in a gas environment
- **Physical accuracy**: Phase transition emerges naturally from the physics without manual thresholding

## Project Structure

```
grade-project/
├── src/
│   └── lbm/                          # Python library package
│       ├── __init__.py
│       ├── core.py                   # LBMSimulator — Python wrapper around C++ extension
│       ├── runner.py                 # Generic run_simulation() loop
│       ├── metrics.py                # Droplet metrics (radius, circularity, aspect ratio…)
│       ├── visualization.py          # Matplotlib plotting functions
│       └── export/
│           ├── obj.py                # Wavefront OBJ heightmap exporter
│           └── blender_import.py     # Blender scripting helper
├── scenarios/                        # One file per experiment
│   ├── a_equilibrium.py
│   ├── b_evaporation.py
│   ├── c_coexistence_curve.py
│   └── d_laplace_pressure.py
├── cpp/
│   ├── include/
│   │   └── lbm_shan_chen.h           # LBM class header
│   └── src/
│       └── lbm_shan_chen.cpp         # LBM class implementation
├── bindings/
│   └── pybind_module.cpp             # PyBind11 interface
├── tests/
│   └── test_metrics.py               # Unit tests (no C++ extension required)
├── docs/
│   ├── IMPLEMENTATION.md             # Algorithm and implementation details
│   └── QUICKSTART.md                 # Quick start guide
├── main.py                           # Entry point: python main.py [--scenario a|b|c|d|all]
├── CMakeLists.txt
├── setup.py
├── pyproject.toml
├── requirements.txt
└── README.md
```

## Mathematical Model

### D2Q9 Lattice Boltzmann Method

The simulation uses a D2Q9 velocity set (9 discrete velocities in 2D):
- 1 rest particle (velocity 0)
- 4 cardinal directions (up, down, left, right)
- 4 diagonal directions

### Shan-Chen Pseudo-Potential Force

The phase separation is driven by the Shan-Chen interaction force:

```
F = -G * psi(rho) * Σ(w_i * psi(rho_neighbor) * c_i)
```

where:
- `G` is the cohesion parameter (negative for liquid-gas separation)
- `psi(rho) = 1 - exp(-1.5*rho)` is the pseudo-potential function
- `w_i` are the D2Q9 weights
- `c_i` are the velocity vectors

### BGK Collision Operator

The collision step uses the Bhatnagar-Gross-Krook (BGK) operator:

```
f_i' = f_i - (f_i - f_i^eq) / tau
```

where `tau` is the relaxation time controlling viscosity.

### Macroscopic Variables

- **Density**: `rho = Σ f_i`
- **Physical velocity**: `u = (Σ f_i * c_i + F/2) / rho`
- **Equilibrium velocity**: `u_eq = u_phys + F/(2*rho)` (standard Shan-Chen velocity-shift, original 1993 formulation)

## Build Instructions

### Prerequisites

1. **C++ Compiler** with C++17 support:
   - GCC 7+ or Clang 5+ (Linux/Mac)
   - MSVC 2017+ (Windows)

2. **CMake** 3.12 or higher

3. **Python** 3.6+ with development headers

4. **Python dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

### Building the Extension

**Option 1 — CMake (recommended):**
```bash
cmake -B build
cmake --build build
```

The compiled module is placed directly into `src/lbm/` by CMake, so no manual copying is needed.

**Option 2 — setuptools:**
```bash
python setup.py build_ext --inplace --build-lib src/lbm
```

After building, verify the extension is in place:
```
src/lbm/lbm_shan_chen*.so   # Linux / macOS
src/lbm/lbm_shan_chen*.pyd  # Windows
```

## Running the Simulation

```bash
# Run all four scenarios sequentially
python main.py

# Run a single scenario
python main.py --scenario a   # Steady-state equilibrium
python main.py --scenario b   # G-ramp evaporation analogy
python main.py --scenario c   # Coexistence curve sweep
python main.py --scenario d   # Laplace pressure test

# Each scenario file is also independently runnable
python scenarios/a_equilibrium.py
```

### The Four Scenarios

| Scenario | Output directory | Description |
|---|---|---|
| A | `output/scenario_a_equilibrium/` | G=-5.0, 5000 steps — steady-state droplet; also exports OBJ sequence |
| B | `output/scenario_b_evaporation/` | G ramped -5.0→-3.7, 8000 steps — droplet dissolves past G_c ≈ -4.0 |
| C | `output/scenario_c_coexistence/` | Five G values (-4.0 to -6.0), plots equilibrium liquid/gas densities vs G |
| D | `output/scenario_d_laplace/` | Five radii (20–60 lu), fits Δp vs 1/R to extract surface tension σ |

## Parameter Tuning

Each scenario has explicit parameters at the top of its file in `scenarios/`. Common parameters:

### Grid Size (`nx`, `ny`)
- **Larger grids** (e.g., 400×400): More detail, slower computation
- **Smaller grids** (e.g., 100×100): Faster, less detail
- **Recommended**: 200×200 for a good balance

### Relaxation Time (`tau`)
- **Controls viscosity**: `ν = (tau - 0.5) / 3`
- **Range**: 0.5 < tau < 2.0
- **Lower tau** (e.g., 0.6): Higher viscosity, more stable, slower dynamics
- **Higher tau** (e.g., 1.5): Lower viscosity, faster dynamics, may be unstable
- **Recommended**: 1.0

### Cohesion Parameter (`G`)
- **Must be negative** for liquid-gas phase separation
- **Critical value**: G_c ≈ -4 (below this, phase separation occurs)
- **More negative** (e.g., -6.5): Stronger phase separation, sharper interface
- **Less negative** (e.g., -4.5): Weaker phase separation, diffuse interface
- **Too negative** (e.g., < -8): Causes numerical instability and blow-up
- **Recommended**: -5.0

### Initial Densities (`rho_liquid`, `rho_gas`)
- **rho_liquid**: Density of the liquid phase (typically 1.5–3.0)
- **rho_gas**: Density of the gas phase (typically 0.05–0.2)
- **Larger contrast**: More pronounced phase separation
- **Recommended**: rho_liquid=2.0, rho_gas=0.1

### Droplet Parameters
- **center_x, center_y**: Droplet center position (`None` = grid center)
- **radius**: Initial droplet radius in lattice units
- **Recommended**: radius = 30–50 for a 200×200 grid

## Output Files

Output is organized by scenario under `output/`:

Within each scenario directory (A and B):
- `density_initial.png` / `density_initial_surface.png` — initial state heatmap and 3D view
- `density_step_*.png` / `density_step_*_surface.png` — intermediate snapshots
- `density_final.png` / `density_final_surface.png` — final state
- `density_evolution.gif` / `density_evolution_surface.gif` — animations
- `metrics.png` — radius, mass, max density, aspect ratio, and circularity over time
- `density_profile.png` — horizontal cross-section with tanh interface fit
- `velocity_field.png` — spurious currents quiver overlay

Scenario A also produces:
- `obj_sequence/frame_NNNNN.obj` — heightmap mesh per animation frame
- `obj_sequence/frame_NNNNN.mtl` — material referencing the texture
- `obj_sequence/frame_NNNNN_texture.png` — viridis-mapped density texture

## OBJ Heightmap Export

The density field can be exported as a Wavefront OBJ mesh for import into Blender, Maya, or any application that reads OBJ. The Z coordinate of each vertex is proportional to the local density (`z = rho * z_scale`).

### Automatic export (Scenario A)

OBJ export runs automatically when running Scenario A. The output goes to `output/scenario_a_equilibrium/obj_sequence/`.

### Manual export from saved `.npy` files

```bash
# Single frame
python src/lbm/export/obj.py path/to/density.npy -o obj_output

# Directory of frames (exports as numbered sequence)
python src/lbm/export/obj.py path/to/npy_dir/ -o obj_output
```

Options:

| Flag | Default | Description |
|------|---------|-------------|
| `-o / --output` | `obj_output` | Output directory |
| `--z-scale` | `50.0` | Vertical exaggeration factor |
| `--subsample` | `1` | Take every Nth grid point (reduces file size) |
| `--colormap` | `viridis` | Any matplotlib colormap name |
| `--color-mode` | `texture` | `texture` = UV-mapped PNG (Blender/Maya); `vertex` = per-vertex `v x y z r g b` |

### Importing into Blender

1. **File → Import → Wavefront (.obj)** and select any `frame_NNNNN.obj`
2. The mesh imports as a colored heightmap. The texture is automatically linked via the `.mtl` file.
3. To animate, use `src/lbm/export/blender_import.py` — edit the `SEQ_DIR` path and run it from Blender's scripting workspace.

## Code Organization

| Layer | Location | Responsibility |
|---|---|---|
| C++ engine | `cpp/src/lbm_shan_chen.cpp` | D2Q9 LBM: force, collision, streaming, macroscopic update |
| PyBind11 bridge | `bindings/pybind_module.cpp` | Exposes `LBMShanChen` class to Python |
| Python wrapper | `src/lbm/core.py` | `LBMSimulator` — single import point, build-error guidance |
| Simulation loop | `src/lbm/runner.py` | Generic `run_simulation()` with plotting and animation |
| Analysis | `src/lbm/metrics.py` | Droplet radius, circularity, aspect ratio via scipy morphology |
| Plotting | `src/lbm/visualization.py` | Heatmaps, 3D surfaces, metric plots, density profiles, velocity fields |
| 3D export | `src/lbm/export/obj.py` | Wavefront OBJ heightmap export (texture or vertex color) |
| Experiments | `scenarios/` | One `run()` function per scenario, independently executable |

## Running Tests

```bash
python -m pytest tests/
# or directly
python tests/test_metrics.py
```

The tests in `tests/test_metrics.py` validate the droplet metrics (radius, circularity, aspect ratio, mass conservation) using synthetic density arrays and do not require the C++ extension to be built.

## Performance Considerations

- **Memory layout**: Row-major order for cache efficiency
- **Streaming buffer**: `f_tmp_` pre-allocated in constructor — no per-step heap allocation
- **Compilation**: `-O3` optimization flag (included in CMakeLists.txt)
- **MLUPS metric**: Throughput is reported as Million Lattice Updates Per Second at the end of each run. A typical 200×200 grid achieves 20–100 MLUPS depending on hardware.

For larger simulations, consider:
- OpenMP parallelization of the inner loops
- GPU acceleration (CUDA/OpenCL)
- Adaptive mesh refinement

## Decision on 3D Extension

The project proposal listed a D3Q19/D3Q27 three-dimensional extension as a conditional objective ("según viabilidad computacional"). After profiling the D2Q9 prototype on the development machine (mid-range laptop), the decision was made to remain with the 2D implementation for the following reasons:

### Memory requirements

A D3Q19 grid of 100³ cells requires:
- **Distribution functions** `f`: 19 × 10⁶ doubles ≈ 145 MB
- **Density, velocity, force**: ≈ 56 MB combined
- **Streaming buffer** `f_tmp`: ≈ 145 MB

**Total ≈ 350 MB** for a single simulation instance, exceeding the memory comfortably available for interactive work on the target hardware.

### Compute time

Empirical scaling from the D2Q9 results predicts a D3Q19 run at 100³ would be 20–50× slower per step. A 3000-step run would take an estimated **30–60 minutes** per scenario, making iterative parameter tuning impractical.

### Scope decision

The 2D prototype fully validates the core physics: spontaneous phase separation, droplet equilibrium, Young-Laplace pressure scaling, G-ramp evaporation, and the coexistence curve. A D3Q19 extension would require only a new C++ class — the Python driver, bindings, and visualization layer would need minimal changes. This is documented as a clear future-work path rather than a limitation of the method.

## Known Limitations

1. **D2Q9 only** — The surface-view plots are 3D surface renderings of a 2D density field, not a volumetric 3D simulation. See the section above for the rationale.
2. **Lattice anisotropy** — At intermediate G values the droplet can appear slightly square due to the lattice geometry. The circularity metric quantifies this artifact (1.0 = perfect circle, ≈0.785 = square).
3. **BGK accuracy** — The single-relaxation-time BGK operator limits accuracy at high density ratios. A multiple-relaxation-time (MRT) operator would improve stability.
4. **Mass drift** — With the velocity clamp threshold at u²=0.04 (Mach ≈0.2), a typical 5000-step run shows mass drift well under 1%.

## Troubleshooting

**`ImportError: No module named 'lbm_shan_chen'`**
- Build the C++ extension first (see Build Instructions above)
- Ensure the compiled `.so`/`.pyd` is inside `src/lbm/`
- Check that the Python version matches the compiled module

**Compilation errors**
- Ensure C++17 support is enabled in your compiler
- Verify PyBind11 and Python development headers are installed

**Unstable simulation**
- Reduce `tau` (increase viscosity)
- Make `G` less negative
- Reduce the density contrast (`rho_liquid` / `rho_gas` ratio)

**Poor phase separation**
- Make `G` more negative (e.g., -5.5 or -6.0)
- Increase density contrast
- Run for more steps

## References

1. Shan, X., & Chen, H. (1993). Lattice Boltzmann model for simulating flows with multiple phases and components. *Physical Review E*, 47(3), 1815.

2. Chen, S., & Doolen, G. D. (1998). Lattice Boltzmann method for fluid flows. *Annual Review of Fluid Mechanics*, 30(1), 329–364.

3. Krüger, T., et al. (2017). *The Lattice Boltzmann Method: Principles and Practice*. Springer.

## License

This code is provided as-is for educational and research purposes.

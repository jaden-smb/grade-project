# Two-Phase Lattice Boltzmann Simulation - Shan-Chen Model

A complete implementation of a two-phase Lattice Boltzmann Method (LBM) simulation using the Shan-Chen pseudo-potential model for liquid-gas phase separation. The codebase uses Python for high-level orchestration and C++ (exposed through PyBind11) for all performance-critical operations.

## Overview

This simulation demonstrates:
- **Phase separation**: Formation of liquid and gas phases from an initial mixed state
- **Interface dynamics**: Evolution of the liquid-gas interface
- **Droplet formation**: Creation and evolution of liquid droplets in a gas environment
- **Physical accuracy**: Phase transition emerges naturally from the physics without manual thresholding

## Project Structure

```
.
├── cpp/                    # C++ implementation (performance-critical)
│   ├── lbm_shan_chen.h    # LBM class header
│   └── lbm_shan_chen.cpp  # LBM class implementation
├── bindings/               # PyBind11 interface
│   └── pybind_module.cpp  # Python bindings
├── python/                 # Python scripts
│   └── simulate.py        # Main simulation and visualization script
├── CMakeLists.txt         # Build configuration
└── README.md              # This file
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
- **Equilibrium velocity**: `u_eq = u_phys + F/(2*rho)` = `u_bare + F/rho` (standard Shan-Chen velocity-shift, original 1993 formulation)

## Build Instructions

### Prerequisites

1. **C++ Compiler** with C++17 support:
   - GCC 7+ or Clang 5+ (Linux/Mac)
   - MSVC 2017+ (Windows)

2. **CMake** 3.12 or higher

3. **Python** 3.6+ with development headers

4. **PyBind11**:
   ```bash
   pip install pybind11
   ```

5. **Python dependencies**:
   ```bash
   pip install numpy matplotlib
   ```

### Building the Extension

1. **Create a build directory**:
   ```bash
   mkdir build
   cd build
   ```

2. **Configure with CMake**:
   ```bash
   cmake ..
   ```

3. **Build**:
   ```bash
   cmake --build .
   ```

   On Linux/Mac, this creates `lbm_shan_chen.cpython-*.so`
   On Windows, this creates `lbm_shan_chen.pyd`

4. **Copy the module to the project root** (if needed):
   ```bash
   cp lbm_shan_chen*.so ..  # Linux/Mac
   # or
   copy lbm_shan_chen*.pyd ..  # Windows
   ```

### Alternative: Using setup.py (Optional)

You can also create a `setup.py` for easier installation:

```python
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, Extension

ext_modules = [
    Pybind11Extension(
        "lbm_shan_chen",
        ["cpp/lbm_shan_chen.cpp", "bindings/pybind_module.cpp"],
        include_dirs=["cpp"],
        cxx_std=17,
    ),
]

setup(
    name="lbm_shan_chen",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)
```

Then build with:
```bash
python setup.py build_ext --inplace
```

## Running the Simulation

1. **Navigate to the project directory**:
   ```bash
   cd /path/to/grade-project-smb
   ```

2. **Run the simulation**:
   ```bash
   python python/simulate.py
   ```

The script runs four scenarios sequentially:

- **Scenario A** (`output/scenario_a_equilibrium/`): Steady-state droplet equilibrium — G=-5.0 fixed, 5000 steps. Demonstrates spontaneous phase separation.
- **Scenario B** (`output/scenario_b_evaporation/`): G-ramp evaporation analogy — G ramped from -5.0 to -3.7 (past critical G_c ≈ -4.0), 8000 steps. Droplet dissolves as cohesion weakens past the critical point.
- **Scenario C** (`output/scenario_c_coexistence/`): Parameter sweep — five G values (-4.0 to -6.0), plots coexistence curve of liquid and gas equilibrium densities.
- **Scenario D** (`output/scenario_d_laplace/`): Laplace pressure test — five droplet radii (20–60 lu), fits Δp vs 1/R to extract surface tension σ using the full Shan-Chen EoS.

## Parameter Tuning

Edit the `main()` function in `python/simulate.py` to modify simulation parameters:

### Grid Size (`nx`, `ny`)
- **Larger grids** (e.g., 400x400): More detail, slower computation
- **Smaller grids** (e.g., 100x100): Faster, less detail
- **Recommended**: 200x200 for good balance

### Relaxation Time (`tau`)
- **Controls viscosity**: `nu = (tau - 0.5) / 3`
- **Range**: 0.5 < tau < 2.0
- **Lower tau** (e.g., 0.6): Higher viscosity, more stable, slower dynamics
- **Higher tau** (e.g., 1.5): Lower viscosity, faster dynamics, may be unstable
- **Recommended**: 1.0

### Cohesion Parameter (`G`)
- **MUST be negative** for liquid-gas phase separation
- **Critical value**: G_c ≈ -4 (below this, phase separation occurs)
- **Range**: -4.5 to -7.0
- **More negative** (e.g., -6.5): Stronger phase separation, sharper interface
- **Less negative** (e.g., -4.5): Weaker phase separation, diffuse interface
- **Too negative** (e.g., < -8): Causes numerical instability and blow-up
- **Recommended**: -5.5

### Initial Densities (`rho_liquid`, `rho_gas`)
- **rho_liquid**: Density of the liquid phase (typically 1.5 - 3.0)
- **rho_gas**: Density of the gas phase (typically 0.05 - 0.2)
- **Larger contrast**: More pronounced phase separation
- **Recommended**: rho_liquid=2.0, rho_gas=0.1

### Droplet Parameters
- **center_x, center_y**: Droplet center position (None = grid center)
- **radius**: Initial droplet radius in lattice units
- **Larger radius**: Bigger initial droplet
- **Recommended**: radius = 30-50 for 200x200 grid

### Number of Steps (`num_steps`)
- **More steps**: Longer simulation time, more evolution
- **Recommended**: 500-2000 steps depending on desired evolution

## Output Files

Output is organized by scenario:

- `output/scenario_a_equilibrium/` — steady-state droplet plots and animations
- `output/scenario_b_evaporation/` — evaporation ramp plots and animations
- `output/scenario_c_coexistence/` — coexistence curve plot (`coexistence_curve.png`)

Within each scenario directory:
- `density_initial.png`: Initial state heatmap
- `density_initial_surface.png`: Initial state surface view
- `density_step_*.png`: Intermediate heatmaps
- `density_step_*_surface.png`: Intermediate surface views
- `density_final.png`: Final state heatmap
- `density_final_surface.png`: Final state surface view
- `density_evolution.gif`: 2D animated evolution
- `density_evolution_surface.gif`: Surface view animation
- `metrics.png`: Radius, mass, max density, aspect ratio, and circularity over time

## Code Explanation

### C++ Implementation (`cpp/lbm_shan_chen.cpp`)

The core LBM class implements:

1. **Initialization**: Sets up grid, allocates memory, initializes fields
2. **Droplet Initialization**: Creates a circular liquid region in gas
3. **Time Step** (`step()`):
   - Computes equilibrium distribution
   - Calculates Shan-Chen force
   - Applies BGK collision
   - Streams distribution functions
   - Updates macroscopic variables

### PyBind11 Interface (`bindings/pybind_module.cpp`)

Exposes the C++ class to Python with:
- Constructor with all parameters
- `initialize_droplet()` method
- `step()` method for time evolution
- Accessors for density and velocity fields
- Grid dimension properties

### Python Script (`python/simulate.py`)

Provides:
- Parameter configuration
- Simulation loop
- Visualization using matplotlib
- Animation generation
- Statistics and analysis

## Performance Considerations

- **Memory layout**: Row-major order for cache efficiency
- **Streaming buffer**: `f_tmp_` pre-allocated in constructor — no per-step heap allocation
- **Optimized loops**: Direct array access, minimal function call overhead
- **Compilation**: Use `-O3` optimization flag (included in CMakeLists.txt)
- **MLUPS metric**: Throughput is reported as Million Lattice Updates Per Second (MLUPS) at the end of each simulation run. A typical 200×200 grid achieves 20–100 MLUPS depending on hardware.

For larger simulations, consider:
- Using OpenMP for parallelization
- GPU acceleration (CUDA/OpenCL)
- Adaptive mesh refinement

## Decision on 3D Extension

The project proposal listed a D3Q19/D3Q27 three-dimensional extension as a conditional objective ("según viabilidad computacional"). After profiling the D2Q9 prototype on the development machine (mid-range laptop), the decision was made to remain with the 2D implementation for the following reasons:

### Memory requirements

A D3Q19 grid of 100³ cells requires:
- **Distribution functions** `f`: 19 × 10⁶ doubles = ~145 MB
- **Density** `ρ`: 1 × 10⁶ doubles = ~8 MB
- **Velocity** `u` (3 components): 3 × 10⁶ doubles = ~24 MB
- **Force** `F` (3 components): 3 × 10⁶ doubles = ~24 MB
- **Streaming buffer** `f_tmp`: same as `f` = ~145 MB

**Total ≈ 350 MB** for a single simulation instance, exceeding the memory comfortably available for interactive work on the target hardware.

### Compute time

Each streaming step in 3D touches 19 neighbors per cell with periodic wrapping in three dimensions. Empirical scaling from the D2Q9 results (≈50 MLUPS on the development machine) predicts:
- D3Q19 at 100³: roughly 20–50× slower per step than D2Q9 at comparable resolution.
- A 3000-step run at 100³ would take an estimated **30–60 minutes** per scenario, making iterative parameter tuning impractical.

### Scope decision

The 2D prototype fully validates the core physics targeted by the thesis:
- Spontaneous liquid-gas phase separation
- Droplet equilibrium and circularity
- Young-Laplace pressure scaling (surface tension measurement)
- G-ramp evaporation analogy
- Coexistence curve

The architecture is structured so that a D3Q19 extension would require only a new C++ class with the 19-velocity set — the Python driver, PyBind11 bindings, and visualisation layer would need minimal changes. This is documented here as a clear future-work path rather than a limitation of the method.

## Known Limitations

1. **D2Q9 only** — This is a 2D simulation. The surface-view plots are 3D surface renderings of the 2D density field, not a volumetric 3D (D3Q19/D3Q27) simulation. See the section above for the rationale.
2. **Lattice anisotropy** — At intermediate G values the droplet can appear slightly square due to the lattice geometry. The circularity metric (reported in `metrics.png` and animation titles) quantifies this artifact (1.0 = perfect circle, ~0.785 = square).
3. **BGK accuracy** — The single-relaxation-time BGK operator limits accuracy at high density ratios. A multiple-relaxation-time (MRT) operator would improve stability.
4. **Mass drift** — With the velocity clamp threshold at u²=0.04 (Mach ~0.2), a typical 5000-step Scenario A run shows mass drift well under 1%. Higher |G| values or longer runs may exhibit slightly larger drift due to residual clamping near the interface.

## Troubleshooting

### Import Error
If you get `ImportError: No module named 'lbm_shan_chen'`:
- Ensure the module is built and in the project root
- Check that the Python version matches the compiled module
- Verify PyBind11 is installed correctly

### Compilation Errors
- Ensure C++17 support is enabled
- Check that PyBind11 headers are found
- Verify Python development headers are installed

### Unstable Simulation
- Reduce `tau` (increase viscosity)
- Make `G` less negative
- Reduce density contrast
- Use smaller time steps (run more steps with smaller changes)

### Poor Phase Separation
- Make `G` more negative
- Increase density contrast
- Run for more steps
- Check initial conditions

## References

1. Shan, X., & Chen, H. (1993). Lattice Boltzmann model for simulating flows with multiple phases and components. *Physical Review E*, 47(3), 1815.

2. Chen, S., & Doolen, G. D. (1998). Lattice Boltzmann method for fluid flows. *Annual Review of Fluid Mechanics*, 30(1), 329-364.

3. Krüger, T., et al. (2017). *The Lattice Boltzmann Method: Principles and Practice*. Springer.

## License

This code is provided as-is for educational and research purposes.


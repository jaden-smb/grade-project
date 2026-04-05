# Implementation Details

## Algorithm Overview

The simulation implements a two-phase LBM using the Shan-Chen pseudo-potential model. The algorithm follows this sequence each time step:

```
1. Compute Shan-Chen force: F = -G * psi(rho) * Σ(w_i * psi(rho_neighbor) * c_i)
2. Compute equilibrium: f_i^eq = w_i * rho * [1 + 3*(c_i·u_eq) + 4.5*(c_i·u_eq)^2 - 1.5*u_eq^2]
   where u_eq = u_phys + F/(2*rho) = u_bare + F/rho
3. Collide: f_i' = f_i - (f_i - f_i^eq) / tau
4. Stream: f_i(x + c_i, t+1) = f_i'(x, t)
5. Update macroscopic: rho = Σ f_i, u = (Σ f_i*c_i + F/2) / rho
```

## Force Implementation

The Shan-Chen force is implemented using the "velocity shift" method:

1. **Force computation**: Based on current density field
   ```cpp
   F = -G * psi(rho) * Σ(w_i * psi(rho_neighbor) * c_i)
   ```

2. **Equilibrium calculation**: Uses the force-shifted velocity
   ```cpp
   u_eq = u + F/(2*rho)   // u already has F/(2*rho) from updateMacroscopic(),
                           // so u_eq = u_bare + F/rho (standard Shan-Chen shift)
   f_i^eq = f_i^eq(u_eq)
   ```

3. **Velocity update**: Accounts for force in momentum
   ```cpp
   u = (Σ f_i*c_i + F/2) / rho
   ```

This approach ensures second-order accuracy in time and maintains stability.

## D2Q9 Velocity Set

The 9 velocity directions are:

```
6  2  5
 \ | /
3--0--1
 / | \
7  4  8
```

Weights:
- Rest particle (0): 4/9
- Cardinal directions (1-4): 1/9
- Diagonal directions (5-8): 1/36

## Boundary Conditions

Currently implemented: **Periodic boundaries**

The streaming step wraps around at grid edges:
- Left edge connects to right edge
- Top edge connects to bottom edge

This is suitable for studying phase separation without wall effects.

## Memory Layout

All fields use row-major order (C-style):
- 1D arrays: `data[j * nx + i]` for 2D position (i, j)
- Distribution functions: `f[(j * nx + i) * Q + q]` for direction q

This layout optimizes cache performance for row-wise access patterns.

## Numerical Stability

The implementation includes several stability measures:

1. **Density threshold**: Avoids division by zero when `rho < 1e-10`
2. **Relaxation time**: Must satisfy `tau > 0.5` for stability
3. **Force magnitude**: Very large |G| can cause instability
4. **Density range**: Extreme densities may cause numerical issues
5. **Velocity clamping counter**: `computeEquilibrium()` clamps velocity when `u² > 0.04` (Mach ~0.2). The clamp count is tracked via `getClampCount()` / `resetClampCount()` and logged in Python every `save_interval` steps. Frequent clamping indicates the simulation is operating near or beyond its stability limit.

## Performance Optimizations

1. **Memory reuse**: Streaming uses temporary buffer to avoid overwriting
2. **Direct array access**: Minimal function call overhead
3. **Compiler optimization**: `-O3` flag enables aggressive optimizations
4. **Cache-friendly layout**: Row-major order for spatial locality

## Extending the Code

### Adding Different Boundary Conditions

Modify the `stream()` function to implement:
- Bounce-back (no-slip walls)
- Free-slip boundaries
- Inflow/outflow boundaries

### Adding Multiple Components

Extend to multi-component Shan-Chen:
- Add separate density fields for each component
- Implement cross-component interaction forces
- Modify equilibrium to account for component interactions

### Parallelization

Potential parallelization strategies:
- **OpenMP**: Parallelize loops over grid points
- **MPI**: Domain decomposition for distributed memory
- **GPU**: CUDA/OpenCL for massive parallelism

### Visualization Enhancements

Add to Python script:
- Velocity field visualization (quiver plots)
- Pressure field computation and display
- Interface tracking and analysis
- Statistical analysis (density distribution, interface width)

## Validation

The implementation can be validated by:

1. **Mass conservation**: Total mass should remain constant. The Python driver reports mass drift percentage at the end of each run (`Mass drift: +X.XXXX%`). After the velocity-clamp threshold was raised to u²=0.04, mass drift in a 5000-step Scenario A run is typically < 1%.
2. **Momentum conservation**: In absence of forces, momentum should be conserved
3. **Phase separation**: For G < 0, system should separate into liquid and gas
4. **Interface properties**: Interface width and surface tension should match theoretical predictions
5. **Coexistence curve sweep**: Scenario C in `simulate.py` sweeps G from -4.5 to -6.5, measuring the final equilibrium liquid (max ρ) and gas (min ρ) densities, and plots the resulting coexistence curve. This provides a quantitative check against Maxwell construction predictions.
6. **Laplace pressure test (Scenario D)**: Δp vs 1/R is fitted linearly to extract surface tension σ. The pressure uses the full Shan-Chen equation of state: `p = ρ·cs² + G·cs²/2·ψ(ρ)²`, where `ψ(ρ) = 1 − exp(−1.5ρ)` and `cs² = 1/3`. This eliminates the large unphysical y-intercept that appears when the ideal-gas EoS (`p = ρ·cs²`) is used instead, and makes the fitted σ physically meaningful.

## Known Limitations

1. **Isothermal**: Temperature is not modeled (single relaxation time)
2. **Single component**: Only one fluid component (two phases of same component)
3. **Periodic boundaries only**: No wall boundaries implemented
4. **2D only**: Extension to 3D requires D3Q19 or D3Q27 velocity set

## References

Key equations and methods are based on:

1. Shan, X., & Chen, H. (1993). Lattice Boltzmann model for simulating flows with multiple phases and components. *Physical Review E*, 47(3), 1815.

2. He, X., Chen, S., & Zhang, R. (1999). A lattice Boltzmann scheme for incompressible multiphase flow and its application in simulation of Rayleigh–Taylor instability. *Journal of Computational Physics*, 152(2), 642-663.

3. Krüger, T., et al. (2017). *The Lattice Boltzmann Method: Principles and Practice*. Springer.


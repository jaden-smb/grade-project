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

The simulation will:
1. Initialize a liquid droplet in a gas environment
2. Run for 1000 time steps
3. Generate plots showing density evolution
4. Create an animated GIF
5. Save output images to the current directory

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


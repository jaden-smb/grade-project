#!/bin/bash
# Build script for LBM Shan-Chen simulation

set -e  # Exit on error

echo "Building LBM Shan-Chen PyBind11 extension..."

# Activate venv if it exists
if [ -d "venv" ]; then
    echo "Activating venv..."
    source venv/bin/activate
elif [ -d ".venv" ]; then
    echo "Activating .venv..."
    source .venv/bin/activate
fi

# Create build directory
mkdir -p build
cd build

# Configure with CMake
echo "Configuring with CMake..."
cmake ..

# Build
echo "Building..."
cmake --build . -j$(nproc)

# Check if module was built (it should already be in project root due to CMake config)
cd ..
if [ -f lbm_shan_chen*.so ] || [ -f lbm_shan_chen*.pyd ]; then
    echo "Build complete! Module is ready in project root."
else
    # Try copying from build directory as fallback
    if [ -f build/lbm_shan_chen*.so ]; then
        cp build/lbm_shan_chen*.so .
        echo "Build complete! Module copied to project root."
    elif [ -f build/lbm_shan_chen*.pyd ]; then
        cp build/lbm_shan_chen*.pyd .
        echo "Build complete! Module copied to project root."
    else
        echo "Warning: Could not find built module."
        exit 1
    fi
fi

echo "Done! You can now run: python python/simulate.py"


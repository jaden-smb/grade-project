cmake_minimum_required(VERSION 3.12)
project(lbm_shan_chen)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Python first (prefer venv Python if available)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# Find PyBind11 using Python to locate it (works with venv)
# This approach finds pybind11 from the Python environment (including venv)
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import pybind11; print(pybind11.get_cmake_dir())"
    OUTPUT_VARIABLE pybind11_CMAKE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE pybind11_FIND_RESULT
    ERROR_QUIET
)

if(pybind11_FIND_RESULT EQUAL 0 AND EXISTS "${pybind11_CMAKE_DIR}/pybind11Config.cmake")
    list(APPEND CMAKE_PREFIX_PATH "${pybind11_CMAKE_DIR}")
    find_package(pybind11 REQUIRED NO_MODULE)
else()
    # Fallback: try standard find_package
    find_package(pybind11 QUIET)
    if(NOT pybind11_FOUND)
        message(FATAL_ERROR 
            "Could not find pybind11. Please install it in your venv with:\n"
            "  source venv/bin/activate\n"
            "  pip install pybind11\n"
            "Or install system package:\n"
            "  sudo apt-get install pybind11-dev  # Ubuntu/Debian\n"
            "  sudo yum install pybind11-devel    # RHEL/CentOS")
    endif()
endif()

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/cpp)

# C++ source files
set(CPP_SOURCES
    cpp/lbm_shan_chen.cpp
)

# PyBind11 module
pybind11_add_module(lbm_shan_chen
    ${CPP_SOURCES}
    bindings/pybind_module.cpp
)

# Compiler-specific options
if(MSVC)
    target_compile_options(lbm_shan_chen PRIVATE /W4)
else()
    target_compile_options(lbm_shan_chen PRIVATE -Wall -Wextra -O3)
endif()

# Set output directory to match Python module location
set_target_properties(lbm_shan_chen PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# Print build information
message(STATUS "Building LBM Shan-Chen module")
message(STATUS "Python executable: ${Python3_EXECUTABLE}")
message(STATUS "Python version: ${Python3_VERSION}")
message(STATUS "PyBind11 version: ${pybind11_VERSION}")


"""
Setup script for building the LBM Shan-Chen PyBind11 extension
Alternative to CMake build method
"""

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, Extension
import pybind11

ext_modules = [
    Pybind11Extension(
        "lbm_shan_chen",
        [
            "cpp/lbm_shan_chen.cpp",
            "bindings/pybind_module.cpp"
        ],
        include_dirs=["cpp"],
        cxx_std=17,
        extra_compile_args=["-O3", "-Wall", "-Wextra"] if pybind11.get_cmake_dir() else [],
    ),
]

setup(
    name="lbm_shan_chen",
    version="1.0.0",
    description="Two-Phase Lattice Boltzmann Method using Shan-Chen Model",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.6",
)


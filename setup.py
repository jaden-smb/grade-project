from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, Extension
import pybind11, sys

ext_modules = [
    Pybind11Extension(
        "lbm_shan_chen",
        [
            "cpp/src/lbm_shan_chen.cpp",
            "bindings/pybind_module.cpp"
        ],
        include_dirs=["cpp/include"],
        cxx_std=17,
        extra_compile_args=[] if sys.platform == "win32" else ["-O3", "-Wall", "-Wextra"],
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


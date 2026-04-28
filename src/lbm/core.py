import sys as _sys
import os as _os
# Ensure the compiled extension is found in the same directory as this file
# (needed when src/ is on sys.path but src/lbm/ is not).
_here = _os.path.dirname(_os.path.abspath(__file__))
if _here not in _sys.path:
    _sys.path.insert(0, _here)

try:
    import lbm_shan_chen as _ext
except ImportError as e:
    raise ImportError(
        "Could not import lbm_shan_chen C++ extension.\n"
        "Build it first with CMake:\n"
        "  cmake -B build && cmake --build build\n"
        "or with setuptools:\n"
        "  pip install pybind11 && python setup.py build_ext --inplace"
    ) from e


class LBMSimulator:

    def __init__(self, nx: int, ny: int, tau: float, G: float,
                 rho_liquid: float, rho_gas: float,
                 eos_type: int = 0, T: float = 0.7 * 0.09433):
        self._sim = _ext.LBMShanChen(nx, ny, tau, G, rho_liquid, rho_gas)
        self.nx = nx
        self.ny = ny
        if eos_type != 0:
            self._sim.set_eos_type(eos_type)
        self._sim.set_temperature(T)

    def initialize_droplet(self, center_x: int, center_y: int, radius: float) -> None:
        self._sim.initialize_droplet(center_x, center_y, radius)

    def step(self) -> None:
        self._sim.step()

    def set_G(self, G: float) -> None:
        self._sim.set_G(G)

    def get_density(self):
        return self._sim.get_density()

    def get_velocity_x(self):
        return self._sim.get_velocity_x()

    def get_velocity_y(self):
        return self._sim.get_velocity_y()

    def get_clamp_count(self) -> int:
        return self._sim.get_clamp_count()

    def reset_clamp_count(self) -> None:
        self._sim.reset_clamp_count()

    def set_eos_type(self, eos_type: int) -> None:
        self._sim.set_eos_type(eos_type)

    def set_temperature(self, T: float) -> None:
        self._sim.set_temperature(T)

    def set_rho0(self, rho0: float) -> None:
        self._sim.set_rho0(rho0)

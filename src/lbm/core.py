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
                 rho_liquid: float, rho_gas: float):
        self._sim = _ext.LBMShanChen(nx, ny, tau, G, rho_liquid, rho_gas)
        self.nx = nx
        self.ny = ny

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

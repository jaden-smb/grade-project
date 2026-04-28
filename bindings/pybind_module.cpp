#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "../cpp/include/lbm_shan_chen.h"

namespace py = pybind11;

PYBIND11_MODULE(lbm_shan_chen, m) {
    m.doc() = "Two-Phase Lattice Boltzmann Method using Shan-Chen Model";
    
    py::class_<LBMShanChen>(m, "LBMShanChen")
        .def(py::init<int, int, double, double, double, double>(),
             "Constructor for LBMShanChen",
             py::arg("nx"), py::arg("ny"), py::arg("tau"), 
             py::arg("G"), py::arg("rho_liquid"), py::arg("rho_gas"))
        
        .def("initialize_droplet", &LBMShanChen::initializeDroplet,
             "Initialize a circular liquid droplet in gas",
             py::arg("center_x"), py::arg("center_y"), py::arg("radius"))
        
        .def("step", &LBMShanChen::step,
             "Perform one LBM time step")
        
        .def("get_density", &LBMShanChen::getDensity,
             "Get density field as a 1D array (row-major order)",
             py::return_value_policy::reference_internal)
        
        .def("get_velocity_x", &LBMShanChen::getVelocityX,
             "Get x-velocity field as a 1D array",
             py::return_value_policy::reference_internal)
        
        .def("get_velocity_y", &LBMShanChen::getVelocityY,
             "Get y-velocity field as a 1D array",
             py::return_value_policy::reference_internal)
        
        .def_property_readonly("nx", &LBMShanChen::getNx,
                               "Grid width")

        .def("set_G", &LBMShanChen::setG,
             "Ramp cohesion to simulate heating / evaporation",
             py::arg("G"))

        .def_property_readonly("G", &LBMShanChen::getG, "Cohesion parameter")
        .def_property_readonly("tau", &LBMShanChen::getTau, "Relaxation time")

        .def("get_clamp_count", &LBMShanChen::getClampCount,
             "Number of cells where velocity was clamped in the last equilibrium computation")
        .def("reset_clamp_count", &LBMShanChen::resetClampCount,
             "Reset the velocity clamp counter")

        .def_property_readonly("ny", &LBMShanChen::getNy,
                               "Grid height")

        // EoS control (default eos_type=0 preserves original SC behaviour)
        .def("set_eos_type", &LBMShanChen::setEosType,
             "Select pseudopotential EoS: 0=original Shan-Chen, 1=Carnahan-Starling",
             py::arg("eos_type"))
        .def_property_readonly("eos_type", &LBMShanChen::getEosType,
                               "Active EoS selector (0=SC, 1=CS)")

        .def("set_temperature", &LBMShanChen::setTemperature,
             "Set temperature for non-ideal EoS (Carnahan-Starling); T < Tc = 0.09433",
             py::arg("T"))
        .def_property_readonly("temperature", &LBMShanChen::getTemperature,
                               "Temperature parameter for non-ideal EoS")

        .def("set_rho0", &LBMShanChen::setRho0,
             "Set rho0 for original SC: psi = 1 - exp(-rho/rho0); default 1/1.5",
             py::arg("rho0"))
        .def_property_readonly("rho0", &LBMShanChen::getRho0,
                               "rho0 parameter for original SC EoS");
}


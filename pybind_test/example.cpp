#include <pybind11/pybind11.h>
#include "lib.hpp"


PYBIND11_MODULE(example, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("test", &out_func, "starting a thread");
    // m.def("mult", &cmult, "A function that multiplies two numbers");
}
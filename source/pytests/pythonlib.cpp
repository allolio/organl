#include "config_ng.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
namespace py = pybind11;
#include "mesh.hpp"
#include "face_eval.hpp"
#include "geocontrol.hpp"
#include "energy.hpp"
#include "nagata.hpp"

PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_MAKE_OPAQUE(Energy);
PYBIND11_MODULE(pynagata, m)
{
  m.doc() = "Python bindings for the Nagata library";
  py::bind_vector<std::vector<double>>(m, "DoubleVector", py::module_local(false));
  py::class_<FaceVertexMesh>(m, "FaceVertexMesh")
    .def(py::init())
    .def( "LoadObj",&FaceVertexMesh::LoadObj )
    .def( "SaveObj",&FaceVertexMesh::SaveObj );

  // py::class_<Gradient>(m, "Gradient");
  
  py::class_<Energy>(m, "Energy");

  py::class_<HelfrichEnergy, Energy>(m, "HelfrichEnergy")
    .def(py::init());

  py::class_<ConstraintHelfrichEnergy, Energy>(m, "ConstraintHelfrichEnergy")
    .def(py::init());

  py::class_<Nagata>(m, "Nagata")
    .def(py::init())
    .def( "FromMesh",&Nagata::FromMesh )
    .def( "ReadFreeze",&Nagata::ReadFreeze )
    .def( "PrepareCoefficientVector",&Nagata::PrepareCoefficientVector )
    .def( "GetCoefficientVector",&Nagata::GetCoefficientVector )
    .def( "SetCoefficientVector",&Nagata::SetCoefficientVector )
    .def( "GetGradientVector", &Nagata::GetGradientVector )
    .def( "TotalArea",&Nagata::TotalArea )
    .def( "TotalEnergy",&Nagata::TotalEnergy )
    .def( "ReadFreeze", &Nagata::ReadFreeze )
    .def( "ReadProperties", &Nagata::ReadProperties )
    .def( "MCStep", &Nagata::MCStep );
}


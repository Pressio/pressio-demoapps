
#ifndef PRESSIODEMOAPPS_PYBINDINGS_MAIN_BINDER_HPP_
#define PRESSIODEMOAPPS_PYBINDINGS_MAIN_BINDER_HPP_

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/iostream.h>
#include <pybind11/stl.h>

#include "types.hpp"
#include "mesh.hpp"
#include "advection.hpp"

PYBIND11_MODULE(MODNAME, mParent)
{

  // --------------------------------
  // cell-centered uniform mesh class
  // --------------------------------
  using ccumesh_t =
    pressiodemoapps::impl::CellCenteredUniformMesh<
      pressiodemoappspy::scalar_t,
      int32_t,
      pressiodemoappspy::py_c_arr,
      pressiodemoappspy::py_c_arr,
      true
    >;

  pybind11::class_<ccumesh_t> meshCl(mParent, "CellCenteredUniformMesh");
  meshCl.def(pybind11::init<std::string>());
  meshCl.def("dimensionality",  &ccumesh_t::dimensionality);
  meshCl.def("stencilMeshSize", &ccumesh_t::stencilMeshSize);
  meshCl.def("sampleMeshSize",  &ccumesh_t::sampleMeshSize);
  meshCl.def("stencilSize",	&ccumesh_t::stencilSize);
  meshCl.def("graph",		&ccumesh_t::graph);
  meshCl.def("dx",		&ccumesh_t::dx);
  meshCl.def("dxInv",		&ccumesh_t::dxInv);
  meshCl.def("dy",		&ccumesh_t::dy);
  meshCl.def("dyInv",		&ccumesh_t::dyInv);
  meshCl.def("viewX",		&ccumesh_t::viewX);
  meshCl.def("viewY",		&ccumesh_t::viewY);

  // -----------------------
  // linead advection 1d
  // -----------------------
  using linadv1d_t =
    pressiodemoapps::ad::impl::LinearAdvT<
      pressiodemoappspy::scalar_t,
      ccumesh_t,
      pressiodemoappspy::py_c_arr,
      pressiodemoappspy::py_c_arr,
      pressiodemoappspy::py_c_arr
    >;

  pybind11::class_<linadv1d_t> ad1dCl(mParent, "LinearAdvection1d");
  ad1dCl.def(pybind11::init<ccumesh_t>());
  ad1dCl.def("totalDofSampleMesh",  &linadv1d_t::totalDofSampleMesh);
  ad1dCl.def("totalDofStencilMesh", &linadv1d_t::totalDofStencilMesh);
  ad1dCl.def("createVelocity",
	     &linadv1d_t::createVelocity,
	     pybind11::return_value_policy::automatic);
  ad1dCl.def("velocity", &linadv1d_t::velocity);
}
#endif

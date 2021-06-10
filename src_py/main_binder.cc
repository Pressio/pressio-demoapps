
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
#include "euler.hpp"

PYBIND11_MODULE(MODNAME, mParent)
{
  using index_t   = int32_t;
  using py_cstyle_arr_int = pybind11::array_t<index_t, pybind11::array::c_style>;
  using py_fstyle_arr_int = pybind11::array_t<index_t, pybind11::array::f_style>;

  // --------------------------------
  // cell-centered uniform mesh class
  // --------------------------------
  using ccumesh_t =
    pressiodemoapps::impl::CellCenteredUniformMesh<
      pressiodemoappspy::scalar_t,
      index_t,
      pressiodemoappspy::py_cstyle_arr_sc,
      py_cstyle_arr_int,
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

  // function that constructs the object directly
  mParent.def("loadCellCenterUniformMesh",
	      &pressiodemoapps::loadCellCenterUniformMesh<ccumesh_t>,
	      pybind11::return_value_policy::take_ownership);

  // ---------------------------------------
  // bind enums to set reconstruction order
  // ---------------------------------------
  pybind11::enum_<pressiodemoapps::reconstructionEnum>(mParent, "reconstructWith")
    .value("firstOrder"    ,
	   pressiodemoapps::reconstructionEnum::firstOrder)
    .value("fifthOrderWeno",
	   pressiodemoapps::reconstructionEnum::fifthOrderWeno);

  // -----------------------
  // linead advection 1d
  // -----------------------
  using linadv1d_t =
    pressiodemoapps::ad::impl::LinearAdvT<
      pressiodemoappspy::scalar_t,
      ccumesh_t,
      pressiodemoappspy::py_cstyle_arr_sc,
      pressiodemoappspy::py_cstyle_arr_sc
    >;

  pybind11::class_<linadv1d_t> ad1dCl(mParent, "LinearAdvection1d");
  ad1dCl.def(pybind11::init<
	     const ccumesh_t &
	     >());

  ad1dCl.def(pybind11::init<
	     const ccumesh_t &,
	     pressiodemoapps::reconstructionEnum
	     >());

  ad1dCl.def("totalDofSampleMesh",  &linadv1d_t::totalDofSampleMesh);
  ad1dCl.def("totalDofStencilMesh", &linadv1d_t::totalDofStencilMesh);
  ad1dCl.def("createVelocity",
	     &linadv1d_t::createVelocity,
	     pybind11::return_value_policy::take_ownership);
  ad1dCl.def("velocity", &linadv1d_t::velocity);

  mParent.def("createPeriodicLinearAdvection1d",
	      &pressiodemoapps::createPeriodicLinearAdvection1dForPy<ccumesh_t, linadv1d_t>,
	      pybind11::return_value_policy::take_ownership);

  // ---------------------------------------
  // bind problem enums
  // ---------------------------------------
  pybind11::enum_<pressiodemoapps::euler1dproblemsEnum>(mParent, "euler1d")
    .value("periodic", pressiodemoapps::euler1dproblemsEnum::periodic)
    .value("sod",      pressiodemoapps::euler1dproblemsEnum::sod)
    .value("lax",      pressiodemoapps::euler1dproblemsEnum::lax);

  pybind11::enum_<pressiodemoapps::euler2dproblemsEnum>(mParent, "euler2d")
    .value("periodic", pressiodemoapps::euler2dproblemsEnum::periodic)
    .value("sedov",    pressiodemoapps::euler2dproblemsEnum::sedov)
    .value("riemann",  pressiodemoapps::euler2dproblemsEnum::riemann);

  // -----------------------
  // Euler 1d
  // -----------------------
  using ee1d_t =
    pressiodemoapps::ee::impl::Euler1dAppT<
      pressiodemoappspy::scalar_t,
      ccumesh_t,
      pressiodemoappspy::py_cstyle_arr_sc,
      pressiodemoappspy::py_cstyle_arr_sc,
      pressiodemoappspy::py_cstyle_arr_sc
    >;

  pybind11::class_<ee1d_t> ee1dCl(mParent, "Euler1dProblem");
  ee1dCl.def(pybind11::init<
	     const ccumesh_t &,
	     const pressiodemoapps::reconstructionEnum,
	     const pressiodemoapps::euler1dproblemsEnum
	     >());

  ee1dCl.def("totalDofSampleMesh",  &ee1d_t::totalDofSampleMesh);
  ee1dCl.def("totalDofStencilMesh", &ee1d_t::totalDofStencilMesh);
  ee1dCl.def("gamma",		     &ee1d_t::gamma);
  ee1dCl.def("initialCondition",
	      &ee1d_t::initialCondition,
	      pybind11::return_value_policy::copy);
  ee1dCl.def("createVelocity",
	      &ee1d_t::createVelocity,
	      pybind11::return_value_policy::automatic);
  ee1dCl.def("velocity", &ee1d_t::velocity);

  mParent.def("createEuler1dProblem",
	      &pressiodemoapps::createEuler1dForPy<ccumesh_t, ee1d_t>,
	      pybind11::return_value_policy::take_ownership);

  // -----------------------
  // Euler 2d
  // -----------------------
  using ee2d_t =
    pressiodemoapps::ee::impl::Euler2dAppT<
      pressiodemoappspy::scalar_t,
    ccumesh_t,
    pressiodemoappspy::py_cstyle_arr_sc, // state type
    pressiodemoappspy::py_cstyle_arr_sc, // velo type
    pressiodemoappspy::py_cstyle_arr_sc  // ghost type
    >;

  pybind11::class_<ee2d_t> ee2dCl(mParent, "Euler2dProblem");
  ee2dCl.def(pybind11::init<
	     const ccumesh_t &,
	     const pressiodemoapps::reconstructionEnum,
	     const pressiodemoapps::euler2dproblemsEnum,
	     const int
	     >());

  ee2dCl.def("totalDofSampleMesh",  &ee2d_t::totalDofSampleMesh);
  ee2dCl.def("totalDofStencilMesh", &ee2d_t::totalDofStencilMesh);
  ee2dCl.def("gamma",		    &ee2d_t::gamma);
  ee2dCl.def("initialCondition",
	      &ee2d_t::initialCondition,
	      pybind11::return_value_policy::copy);
  ee2dCl.def("createVelocity",
	      &ee2d_t::createVelocity,
	      pybind11::return_value_policy::automatic);
  ee2dCl.def("velocity", &ee2d_t::velocity);

  mParent.def("createEuler2dProblem",
	      &pressiodemoapps::createEuler2dForPy<ccumesh_t, ee2d_t>,
	      pybind11::return_value_policy::take_ownership);

}
#endif

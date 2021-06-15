
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
#include "euler1d.hpp"
#include "euler2d.hpp"

// helper function to bind API methods
// since they are the same for any application problem
namespace pdabindimpl{
template<class cpp_t, class py_t>
void bindCommonApiMethods(py_t & appObj)
{
  appObj.def("totalDofSampleMesh",  &cpp_t::totalDofSampleMesh);
  appObj.def("totalDofStencilMesh", &cpp_t::totalDofStencilMesh);

  appObj.def("initialCondition",
	     &cpp_t::initialCondition,
	     pybind11::return_value_policy::copy);

  appObj.def("createVelocity",
	     &cpp_t::createVelocity,
	     pybind11::return_value_policy::automatic);

  appObj.def("velocity", &cpp_t::velocity);
}
};

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
      // ordinal type
      index_t,
      // coordinate storage type
      pressiodemoappspy::py_cstyle_arr_sc,
      // graph type
      py_cstyle_arr_int,
      // true because we are doing bindings
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
  meshCl.def("dy",		&ccumesh_t::dy);
  meshCl.def("dxInv",		&ccumesh_t::dxInv);
  meshCl.def("dyInv",		&ccumesh_t::dyInv);
  meshCl.def("viewX",		&ccumesh_t::viewX);
  meshCl.def("viewY",		&ccumesh_t::viewY);

  // function that constructs the object directly
  mParent.def("loadCellCenterUniformMesh",
	      &pressiodemoapps::loadCellCenterUniformMesh<ccumesh_t>,
	      pybind11::return_value_policy::take_ownership);

  // ---------------------------------------
  // bind enums
  // ---------------------------------------
  pybind11::enum_<pressiodemoapps::ReconstructionType>(mParent, "ReconstructionType")
    .value("firstOrder",     pressiodemoapps::ReconstructionType::firstOrder)
    .value("fifthOrderWeno", pressiodemoapps::ReconstructionType::fifthOrderWeno);

  pybind11::enum_<pressiodemoapps::FluxType>(mParent, "FluxType")
    .value("Rusanov",     pressiodemoapps::FluxType::Rusanov);

  pybind11::enum_<pressiodemoapps::Euler1d>(mParent, "Euler1d")
    .value("PeriodicSmooth", pressiodemoapps::Euler1d::PeriodicSmooth)
    .value("Sod",	     pressiodemoapps::Euler1d::Sod)
    .value("Lax",	     pressiodemoapps::Euler1d::Lax);

  pybind11::enum_<pressiodemoapps::Euler2d>(mParent, "Euler2d")
    .value("PeriodicSmooth", pressiodemoapps::Euler2d::PeriodicSmooth)
    .value("SedovFull",	     pressiodemoapps::Euler2d::SedovFull)
    .value("SedovSymmtry",   pressiodemoapps::Euler2d::SedovSymmetry)
    .value("Riemann",	     pressiodemoapps::Euler2d::Riemann);


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

  pybind11::class_<linadv1d_t> adDiff1dProb(mParent, "LinearAdvection1dProblem");
  // adDiff1dProb.def(pybind11::init<
  // 		   const ccumesh_t &,
  // 		   pressiodemoapps::ReconstructionType
  // 		   >());

  pdabindimpl::bindCommonApiMethods<linadv1d_t>(adDiff1dProb);

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
  pdabindimpl::bindCommonApiMethods<ee1d_t>(ee1dCl);
  ee1dCl.def("gamma",		     &ee1d_t::gamma);

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
  pdabindimpl::bindCommonApiMethods<ee2d_t>(ee2dCl);
  ee2dCl.def("gamma",		    &ee2d_t::gamma);

  // -----------------------------------------------
  // functions to create problems.
  // add more as the c++ grows
  // -----------------------------------------------
  mParent.def("createPeriodicLinearAdvection1d",
	      &pressiodemoapps::createPeriodicLinearAdvection1dForPy<ccumesh_t, linadv1d_t>,
	      pybind11::return_value_policy::take_ownership);

  mParent.def("createProblem",
	      &pressiodemoapps::createEuler1dForPyA<ccumesh_t, ee1d_t>,
	      pybind11::return_value_policy::take_ownership);

  mParent.def("createProblem",
	      &pressiodemoapps::createEuler1dForPyB<ccumesh_t, ee1d_t>,
	      pybind11::return_value_policy::take_ownership);

  mParent.def("createProblem",
	      &pressiodemoapps::createEuler2dForPyA<ccumesh_t, ee2d_t>,
	      pybind11::return_value_policy::take_ownership);

  mParent.def("createProblem",
	      &pressiodemoapps::createEuler2dForPyB<ccumesh_t, ee2d_t>,
	      pybind11::return_value_policy::take_ownership);

  mParent.def("createProblem",
	      &pressiodemoapps::createEuler2dForPyC<ccumesh_t, ee2d_t>,
	      pybind11::return_value_policy::take_ownership);

}
#endif

/*
//@HEADER
// ************************************************************************
//
// main_binder.cc
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef PRESSIODEMOAPPS_PYBINDINGS_MAIN_BINDER_HPP_
#define PRESSIODEMOAPPS_PYBINDINGS_MAIN_BINDER_HPP_

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/iostream.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "types.hpp"
#include "pressiodemoapps/mesh.hpp"
// 1d
#include "pressiodemoapps/advection1d.hpp"
#include "pressiodemoapps/diffusion_reaction1d.hpp"
#include "pressiodemoapps/euler1d.hpp"
// 2d
#include "pressiodemoapps/advection_diffusion2d.hpp"
#include "pressiodemoapps/diffusion_reaction2d.hpp"
#include "pressiodemoapps/euler2d.hpp"
#include "pressiodemoapps/swe2d.hpp"
// 3d
#include "pressiodemoapps/euler3d.hpp"

namespace pressiodemoappspy{
namespace impl{

template<class p_t, class T>
void bindCommonApiMethods(T & appObj)
{
  appObj.def("totalDofSampleMesh",
	     &p_t::totalDofSampleMesh);
  appObj.def("totalDofStencilMesh",
	     &p_t::totalDofStencilMesh);

  appObj.def("initialCondition",
	     &p_t::initialCondition,
	     pybind11::return_value_policy::copy);

  appObj.def("createRightHandSide",
	     &p_t::createVelocity,
	     pybind11::return_value_policy::automatic);

  appObj.def("rightHandSide",
	     &p_t::velocity);

  appObj.def("createApplyJacobianResult",
	     &p_t::createApplyJacobianResultRank1);
  appObj.def("createApplyJacobianResult",
	     &p_t::createApplyJacobianResultRank2_ll);
  appObj.def("createApplyJacobianResult",
	     &p_t::createApplyJacobianResultRank2_lr);
  appObj.def("applyJacobian",
	     &p_t::applyJacobianRank1);
  appObj.def("applyJacobian",
	     &p_t::applyJacobianRank2_ll);
  appObj.def("applyJacobian",
	     &p_t::applyJacobianRank2_lr);
}

// ---------------------------
void bindSchemeEnums(pybind11::module & mParent)
{
  namespace pda = pressiodemoapps;

  // reconstructions
  pybind11::enum_<pda::InviscidFluxReconstruction>(mParent, "InviscidFluxReconstruction")
    .value("FirstOrder",
	   pda::InviscidFluxReconstruction::FirstOrder)
    .value("Weno3",
	   pda::InviscidFluxReconstruction::Weno3)
    .value("Weno5",
	   pda::InviscidFluxReconstruction::Weno5);

  pybind11::enum_<pda::ViscousFluxReconstruction>(mParent, "ViscousFluxReconstruction")
    .value("FirstOrder",
	   pda::ViscousFluxReconstruction::FirstOrder);

  // flux schemes
  pybind11::enum_<pda::InviscidFluxScheme>(mParent, "InviscidFluxScheme")
    .value("Rusanov",
	   pda::InviscidFluxScheme::Rusanov);

  pybind11::enum_<pda::ViscousFluxScheme>(mParent, "ViscousFluxScheme")
    .value("Central",
	   pda::ViscousFluxScheme::Central);
}

void bind1dProblemEnums(pybind11::module & mParent)
{
  namespace pda = pressiodemoapps;

  pybind11::enum_<pda::Advection1d>(mParent, "Advection1d")
    .value("PeriodicLinear",
	   pda::Advection1d::PeriodicLinear);

  pybind11::enum_<pda::DiffusionReaction1d>(mParent, "DiffusionReaction1d")
    .value("ProblemA",
	   pda::DiffusionReaction1d::ProblemA);

  pybind11::enum_<pda::Euler1d>(mParent, "Euler1d")
    .value("PeriodicSmooth",
	   pda::Euler1d::PeriodicSmooth)
    .value("Sod",
	   pda::Euler1d::Sod)
    .value("Lax",
	   pda::Euler1d::Lax)
    .value("ShuOsher",
	   pda::Euler1d::ShuOsher);
}

void bind2dProblemEnums(pybind11::module & mParent)
{
  namespace pda = pressiodemoapps;

  pybind11::enum_<pda::AdvectionDiffusion2d>(mParent, "AdvectionDiffusion2d")
    .value("BurgersPeriodic",
	   pda::AdvectionDiffusion2d::BurgersPeriodic);

  pybind11::enum_<pda::DiffusionReaction2d>(mParent, "DiffusionReaction2d")
    .value("ProblemA",
	   pda::DiffusionReaction2d::ProblemA)
    .value("GrayScott",
	   pda::DiffusionReaction2d::GrayScott);

  pybind11::enum_<pda::Euler2d>(mParent, "Euler2d")
    .value("PeriodicSmooth",
	   pda::Euler2d::PeriodicSmooth)
    .value("KelvinHelmholtz",
	   pda::Euler2d::KelvinHelmholtz)
    .value("SedovFull",
	   pda::Euler2d::SedovFull)
    .value("SedovSymmetry",
	   pda::Euler2d::SedovSymmetry)
    .value("Riemann",
	   pda::Euler2d::Riemann)
    .value("NormalShock",
	   pda::Euler2d::NormalShock)
    .value("DoubleMachReflection",
	   pda::Euler2d::DoubleMachReflection);

  pybind11::enum_<pda::Swe2d>(mParent, "Swe2d")
    .value("SlipWall", pda::Swe2d::SlipWall);
}

void bind3dProblemEnums(pybind11::module & mParent)
{
  namespace pda = pressiodemoapps;
  pybind11::enum_<pda::Euler3d>(mParent, "Euler3d")
    .value("PeriodicSmooth", pda::Euler3d::PeriodicSmooth)
    .value("SedovSymmetry",  pda::Euler3d::SedovSymmetry);
}

// ---------------------------
struct CcuMeshBinder
{
  using mesh_t =
    pressiodemoapps::impl::CellCenteredUniformMesh<
      pressiodemoappspy::scalar_t,
      // ordinal type
      pressiodemoappspy::ordinal_t,
      // coordinate storage type
      pressiodemoappspy::py_cstyle_arr_sc,
      // graph type
      pressiodemoappspy::py_cstyle_arr_int,
      // true because we are doing bindings
      true
    >;

  void operator()(pybind11::module & mParent)
  {
    namespace pda = pressiodemoapps;
    pybind11::class_<mesh_t> meshClass(mParent, "CellCenteredUniformMesh");
    meshClass.def(pybind11::init<std::string>());
    meshClass.def("dimensionality",  &mesh_t::dimensionality);
    meshClass.def("stencilMeshSize", &mesh_t::stencilMeshSize);
    meshClass.def("sampleMeshSize",  &mesh_t::sampleMeshSize);
    meshClass.def("stencilSize",     &mesh_t::stencilSize);
    meshClass.def("graph",	     &mesh_t::graph);
    meshClass.def("dx",		     &mesh_t::dx);
    meshClass.def("dy",		     &mesh_t::dy);
    meshClass.def("dz",		     &mesh_t::dz);
    meshClass.def("dxInv",	     &mesh_t::dxInv);
    meshClass.def("dyInv",	     &mesh_t::dyInv);
    meshClass.def("dzInv",	     &mesh_t::dzInv);
    meshClass.def("viewX",	     &mesh_t::viewX);
    meshClass.def("viewY",	     &mesh_t::viewY);
    meshClass.def("viewZ",	     &mesh_t::viewZ);

    // function that constructs the object directly
    mParent.def("load_cellcentered_uniform_mesh",
		&pda::loadCellCenterUniformMesh<mesh_t>,
		pybind11::return_value_policy::take_ownership);
  }
};

// -----------------------
// advection 1d
// -----------------------
template<class MeshType>
void bindAdvection1d(pybind11::module & mParent)
{
  namespace pda = pressiodemoapps;

  using py_problem_type = pda::PublicProblemMixinPy<
    pda::impladvection1d::EigenApp<MeshType>>;

  pybind11::class_<py_problem_type> prob(mParent, "Advection1dProblem");
  pressiodemoappspy::impl::bindCommonApiMethods<py_problem_type>(prob);

  mParent.def("create_problem",
	      &pda::create_linear_advection_1d_problem_default_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert());
  mParent.def("create_linear_advection_1d_problem",
	      &pda::create_linear_advection_1d_problem_ov1_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert());
  mParent.def("create_linear_advection_1d_problem",
	      &pda::create_linear_advection_1d_problem_ov2_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg("ic") = 1);
}

// -----------------------
// diffusion-reaction 1d
// -----------------------
template<class MeshType>
void bindDiffusionReaction1d(pybind11::module & mParent)
{
  namespace pda = pressiodemoapps;

  using py_problem_type = pda::PublicProblemMixinPy<
    pda::impldiffusionreaction1d::EigenApp<MeshType>>;

  pybind11::class_<py_problem_type> prob(mParent, "DiffusionReaction1dProblem");
  pressiodemoappspy::impl::bindCommonApiMethods<py_problem_type>(prob);

  mParent.def("create_problem",
	      &pda::create_diffusion_reaction_1d_problem_default_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert());
  mParent.def("create_diffusion_reaction_1d_problem_A",
	      &pda::create_diffusion_reaction_1d_problem_A_ov1_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert());
  mParent.def("create_diffusion_reaction_1d_problem_A",
	      &pda::create_diffusion_reaction_1d_problem_A_ov2_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert());
}

// -----------------------
// euler 1d
// -----------------------
template<class MeshType>
void bindEuler1d(pybind11::module & mParent)
{
  namespace pda = pressiodemoapps;

  using py_problem_type = pda::PublicProblemMixinPy<
    pda::impleuler1d::EigenApp<MeshType>>;

  pybind11::class_<py_problem_type> prob(mParent, "Euler1dProblem");
  pressiodemoappspy::impl::bindCommonApiMethods<py_problem_type>(prob);
  prob.def("gamma", &py_problem_type::gamma);

  mParent.def("create_problem",
	      &pda::create_euler_1d_py_problem_default<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert());
}

// -----------------------
// advection-diffusion 2d
// -----------------------
template<class MeshType>
void bindAdvectionDiffusion2d(pybind11::module & mParent)
{
  namespace pda = pressiodemoapps;

  using py_problem_type = pda::PublicProblemMixinPy<
    pda::impladvdiff2d::EigenApp<MeshType>>;

  pybind11::class_<py_problem_type> prob(mParent, "AdvectionDiffusion2dProblem");
  pressiodemoappspy::impl::bindCommonApiMethods<py_problem_type>(prob);

  mParent.def("create_problem",
	      &pda::create_advecdiffusion2d_problem_default_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert());
  mParent.def("create_periodic_burgers_2d_problem",
	      &pda::create_periodic_burgers_2d_problem_ov1_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(), pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(), pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(), pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(), pybind11::arg().noconvert());
}

// -----------------------
// diffusion-reaction 2d
// -----------------------
template<class MeshType>
void bindDiffusionReaction2d(pybind11::module & mParent)
{
  namespace pda = pressiodemoapps;

  using py_problem_type = pda::PublicProblemMixinPy<
    pda::impldiffusionreaction2d::EigenApp<MeshType>>;

  pybind11::class_<py_problem_type> prob(mParent, "DiffusionReaction2dProblem");
  pressiodemoappspy::impl::bindCommonApiMethods<py_problem_type>(prob);

  mParent.def("create_problem",
	      &pda::create_diffreac2d_problem_default_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert());
  mParent.def("create_problem",
	      &pda::create_diffreac2d_problem_ov1_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert());

  mParent.def("create_diffusion_reaction_2d_problem_A",
	      &pda::create_diffreac2d_problem_A_ov1_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(), pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(), pybind11::arg().noconvert());
  mParent.def("create_diffusion_reaction_2d_problem_A",
	      &pda::create_diffreac2d_problem_A_ov2_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(), pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(), pybind11::arg().noconvert(),
	      pybind11::arg().noconvert());

  mParent.def("create_gray_scott_2d_problem",
	      &pda::create_diffreac2d_problem_grayscott_ov1_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(), pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(), pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(), pybind11::arg().noconvert());
}

// -----------------------
// euler 2d
// -----------------------
template<class MeshType>
void bindEuler2d(pybind11::module & mParent)
{
  namespace pda = pressiodemoapps;

  using py_problem_type = pda::PublicProblemMixinPy<
    pda::impleuler2d::EigenApp<MeshType>>;

  pybind11::class_<py_problem_type> prob(mParent, "Euler2dProblem");
  pressiodemoappspy::impl::bindCommonApiMethods<py_problem_type>(prob);
  prob.def("gamma", &py_problem_type::gamma);

  mParent.def("create_problem",
	      &pda::create_euler_2d_problem_default_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert());
  mParent.def("create_problem",
	      &pda::create_euler_2d_problem_ov1_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert());
  mParent.def("create_cross_shock_problem",
	      &pda::create_euler_2d_cross_shock_problem_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert());
}

// -----------------------
// swe 2d
// -----------------------
template<class MeshType>
void bindSwe2d(pybind11::module & mParent)
{
  namespace pda = pressiodemoapps;

  using py_problem_type = pda::PublicProblemMixinPy<
    pda::implswe2d::EigenApp<MeshType>>;

  pybind11::class_<py_problem_type> prob(mParent, "Swe2dProblem");
  pressiodemoappspy::impl::bindCommonApiMethods<py_problem_type>(prob);
  prob.def("coriolis", &py_problem_type::coriolis);
  prob.def("gravity",  &py_problem_type::gravity);

  mParent.def("create_problem",
	      &pda::create_swe2d_problem_default_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert());
  mParent.def("create_slip_wall_swe_2d_problem",
	      &pda::create_slip_wall_swe2d_problem_ov1_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(), pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(), pybind11::arg().noconvert(),
	      pybind11::arg().noconvert());
}

// -----------------------
// euler 3d
// -----------------------
template<class MeshType>
void bindEuler3d(pybind11::module & mParent)
{
  namespace pda = pressiodemoapps;

  using py_problem_type = pda::PublicProblemMixinPy<
    pda::impleuler3d::EigenApp<MeshType>>;

  pybind11::class_<py_problem_type> prob(mParent, "Euler3dProblem");
  pressiodemoappspy::impl::bindCommonApiMethods<py_problem_type>(prob);
  prob.def("gamma", &py_problem_type::gamma);

  mParent.def("create_problem",
	      &pda::create_euler_3d_problem_default_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert());
  mParent.def("create_problem",
	      &pda::create_euler_3d_problem_ov1_for_py<MeshType, py_problem_type>,
	      pybind11::return_value_policy::take_ownership,
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert(),
	      pybind11::arg().noconvert());
}

}}//end namespace pressiodemoappspy::impl

//=======================================
PYBIND11_MODULE(MODNAME, mTopLevel)
//=======================================
{
  namespace pdapyimpl = pressiodemoappspy::impl;

  pdapyimpl::bindSchemeEnums(mTopLevel);
  pdapyimpl::bind1dProblemEnums(mTopLevel);
  pdapyimpl::bind2dProblemEnums(mTopLevel);
  pdapyimpl::bind3dProblemEnums(mTopLevel);

  // ---------------------------------
  // cell-centered uniform (ccu) mesh
  // ---------------------------------
  using mesh_binder_t = pdapyimpl::CcuMeshBinder;
  using ccumesh_t = typename mesh_binder_t::mesh_t;
  mesh_binder_t ccuMeshB;
  ccuMeshB(mTopLevel);

  pdapyimpl::bindAdvection1d<ccumesh_t>(mTopLevel);
  pdapyimpl::bindDiffusionReaction1d<ccumesh_t>(mTopLevel);
  pdapyimpl::bindEuler1d<ccumesh_t>(mTopLevel);
  pdapyimpl::bindAdvectionDiffusion2d<ccumesh_t>(mTopLevel);
  pdapyimpl::bindDiffusionReaction2d<ccumesh_t>(mTopLevel);
  pdapyimpl::bindEuler2d<ccumesh_t>(mTopLevel);
  pdapyimpl::bindSwe2d<ccumesh_t>(mTopLevel);
  pdapyimpl::bindEuler3d<ccumesh_t>(mTopLevel);

}
#endif


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
#include "pressiodemoapps/advection.hpp"
#include "pressiodemoapps/diffusion_reaction.hpp"
#include "pressiodemoapps/euler1d.hpp"
#include "pressiodemoapps/euler2d.hpp"
#include "pressiodemoapps/swe2d.hpp"
#include "pressiodemoapps/euler3d.hpp"

namespace pressiodemoappspy{ namespace impl{

template<class p_t, class T>
void bindCommonApiMethods(T & appObj)
{
  appObj.def("totalDofSampleMesh",  &p_t::totalDofSampleMesh);
  appObj.def("totalDofStencilMesh", &p_t::totalDofStencilMesh);

  appObj.def("initialCondition",
	     &p_t::initialCondition,
	     pybind11::return_value_policy::copy);

  appObj.def("createVelocity",
	     &p_t::createVelocity,
	     pybind11::return_value_policy::automatic);

  appObj.def("velocity", &p_t::velocity);

  appObj.def("createApplyJacobianResult", &p_t::createApplyJacobianResultRank1);
  appObj.def("createApplyJacobianResult", &p_t::createApplyJacobianResultRank2_ll);
  appObj.def("createApplyJacobianResult", &p_t::createApplyJacobianResultRank2_lr);
  appObj.def("applyJacobian", &p_t::applyJacobianRank1);
  appObj.def("applyJacobian", &p_t::applyJacobianRank2_ll);
  appObj.def("applyJacobian", &p_t::applyJacobianRank2_lr);
}

// ---------------------------
void bindSchemeEnums(pybind11::module & mParent)
{
  namespace pda = pressiodemoapps;

  pybind11::enum_<pda::InviscidFluxReconstruction>(mParent, "InviscidFluxReconstruction")
    .value("FirstOrder",     pda::InviscidFluxReconstruction::FirstOrder)
    .value("Weno3",	     pda::InviscidFluxReconstruction::Weno3)
    .value("Weno5",	     pda::InviscidFluxReconstruction::Weno5);

  pybind11::enum_<pda::ViscousFluxReconstruction>(mParent, "ViscousFluxReconstruction")
    .value("FirstOrder",     pda::ViscousFluxReconstruction::FirstOrder);

  pybind11::enum_<pda::InviscidFluxScheme>(mParent, "InviscidFluxScheme")
    .value("Rusanov",	     pda::InviscidFluxScheme::Rusanov);

  pybind11::enum_<pda::ViscousFluxScheme>(mParent, "ViscousFluxScheme")
    .value("Central",	     pda::ViscousFluxScheme::Central);
}

void bindProblemEnums(pybind11::module & mParent)
{
  namespace pda = pressiodemoapps;

  pybind11::enum_<pda::Advection1d>(mParent, "Advection1d")
    .value("PeriodicLinear", pda::Advection1d::PeriodicLinear);

  pybind11::enum_<pda::DiffusionReaction1d>(mParent, "DiffusionReaction1d")
    .value("ProblemA", pda::DiffusionReaction1d::ProblemA);

  pybind11::enum_<pda::DiffusionReaction2d>(mParent, "DiffusionReaction2d")
    .value("ProblemA", pda::DiffusionReaction2d::ProblemA)
    .value("GrayScott", pda::DiffusionReaction2d::GrayScott);

  pybind11::enum_<pda::Euler1d>(mParent, "Euler1d")
    .value("PeriodicSmooth", pda::Euler1d::PeriodicSmooth)
    .value("Sod",	     pda::Euler1d::Sod)
    .value("Lax",	     pda::Euler1d::Lax)
    .value("ShuOsher",	     pda::Euler1d::ShuOsher);

  pybind11::enum_<pda::Euler2d>(mParent, "Euler2d")
    .value("PeriodicSmooth", pda::Euler2d::PeriodicSmooth)
    .value("KelvinHelmholtz", pda::Euler2d::KelvinHelmholtz)
    .value("SedovFull",	     pda::Euler2d::SedovFull)
    .value("SedovSymmetry",  pda::Euler2d::SedovSymmetry)
    .value("Riemann",	     pda::Euler2d::Riemann)
    .value("NormalShock",    pda::Euler2d::NormalShock)
    .value("RayleighTaylor", pda::Euler2d::RayleighTaylor)
    .value("DoubleMachReflection", pda::Euler2d::DoubleMachReflection);

  pybind11::enum_<pda::Swe2d>(mParent, "Swe2d")
    .value("SlipWall", pda::Swe2d::SlipWall);

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
    meshClass.def("dxInv",	     &mesh_t::dxInv);
    meshClass.def("dyInv",	     &mesh_t::dyInv);
    meshClass.def("viewX",	     &mesh_t::viewX);
    meshClass.def("viewY",	     &mesh_t::viewY);

    // function that constructs the object directly
    mParent.def("load_cellcentered_uniform_mesh",
		&pda::loadCellCenterUniformMesh<mesh_t>,
		pybind11::return_value_policy::take_ownership);
  }
};
}}//end namespace pressiodemoappspy::impl

//=======================================
PYBIND11_MODULE(MODNAME, mTopLevel)
//=======================================
{
  namespace pda = pressiodemoapps;

  pressiodemoappspy::impl::bindSchemeEnums(mTopLevel);
  pressiodemoappspy::impl::bindProblemEnums(mTopLevel);

  // ---------------------------------
  // cell-centered uniform (ccu) mesh
  // ---------------------------------
  using mesh_binder_t = pressiodemoappspy::impl::CcuMeshBinder;
  using ccumesh_t = typename mesh_binder_t::mesh_t;
  mesh_binder_t ccuMeshB;
  ccuMeshB(mTopLevel);

  // -----------------------
  // advection 1d
  // -----------------------
  using ad1d_t = decltype
    (
     pda::impladv::create_problem_for_py
     (
      std::declval<const ccumesh_t &>(),
      std::declval<pda::Advection1d>(),
      std::declval<pda::InviscidFluxReconstruction>()
      )
     );

  pybind11::class_<ad1d_t> adv1dProb(mTopLevel, "Advection1dProblem");
  pressiodemoappspy::impl::bindCommonApiMethods<ad1d_t>(adv1dProb);

  mTopLevel.def("create_problem",
		&pda::impladv::create_problem_for_py<ccumesh_t>,
		pybind11::return_value_policy::take_ownership);

  // -----------------------
  // diffusion-reaction 1d
  // -----------------------
  using diffreac_impl_1d_t = pda::impldiffreac::EigenDiffReac1dApp<ccumesh_t>;

  using diffreac_1d_t = pda::PublicProblemMixinPy<diffreac_impl_1d_t>;
  pybind11::class_<diffreac_1d_t> diffReac1dProb(mTopLevel, "DiffusionReaction1dProblem");
  pressiodemoappspy::impl::bindCommonApiMethods<diffreac_1d_t>(diffReac1dProb);

  mTopLevel.def("create_problem",
		&pda::impldiffreac::create_problem_for_pyA<
		ccumesh_t, diffreac_1d_t, pda::DiffusionReaction1d>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("create_problem",
		&pda::impldiffreac::create_problem_for_pyB<
		ccumesh_t, diffreac_1d_t, pda::DiffusionReaction1d>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("create_problem",
		&pda::impldiffreac::create_problem_for_pyC1d<ccumesh_t, diffreac_1d_t>,
		pybind11::return_value_policy::take_ownership);

  // -----------------------
  // Euler 1d
  // -----------------------
  using euler_impl_1d_t = pda::implee1d::EigenEuler1dApp<ccumesh_t>;

  using euler_1d_t = pressiodemoapps::PublicProblemMixinPy<euler_impl_1d_t>;
  pybind11::class_<euler_1d_t> ee1dProb(mTopLevel, "Euler1dProblem");
  pressiodemoappspy::impl::bindCommonApiMethods<euler_1d_t>(ee1dProb);
  ee1dProb.def("gamma", &euler_1d_t::gamma);

  mTopLevel.def("create_problem",
		&pda::implee1d::create_problem_for_pyA<ccumesh_t, euler_1d_t>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("create_problem",
		&pda::implee1d::create_problem_for_pyB<ccumesh_t, euler_1d_t>,
		pybind11::return_value_policy::take_ownership);

  // -----------------------
  // diffusion-reaction 2d
  // -----------------------
  using diffreac_impl_2d_t = pda::impldiffreac::EigenDiffReac2dApp<ccumesh_t>;

  using diffreac_2d_t = pda::PublicProblemMixinPy<diffreac_impl_2d_t>;
  pybind11::class_<diffreac_2d_t> diffReac2dProb(mTopLevel, "DiffusionReaction2dProblem");
  pressiodemoappspy::impl::bindCommonApiMethods<diffreac_2d_t>(diffReac2dProb);

  mTopLevel.def("create_problem",
		&pda::impldiffreac::create_problem_for_pyA<
		ccumesh_t, diffreac_2d_t, pda::DiffusionReaction2d>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("create_problem",
		&pda::impldiffreac::create_problem_for_pyB<
		ccumesh_t, diffreac_2d_t, pda::DiffusionReaction2d>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("create_problem",
		&pda::impldiffreac::create_problem_for_pyC2d<ccumesh_t, diffreac_2d_t>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("create_problem",
		&pda::impldiffreac::create_problem_for_pyD<
		ccumesh_t, diffreac_2d_t, pda::DiffusionReaction2d>,
		pybind11::return_value_policy::take_ownership);

  // -----------------------
  // Euler 2d
  // -----------------------
  using euler_impl_2d_t = pda::ee::impl::EigenEuler2dApp<ccumesh_t>;

  using euler_2d_t = pressiodemoapps::PublicProblemMixinPy<euler_impl_2d_t>;
  pybind11::class_<euler_2d_t> ee2dProb(mTopLevel, "Euler2dProblem");
  pressiodemoappspy::impl::bindCommonApiMethods<euler_2d_t>(ee2dProb);
  ee2dProb.def("gamma", &euler_2d_t::gamma);

  mTopLevel.def("create_problem",
		&pda::implee2d::create_problem_for_pyA<ccumesh_t, euler_2d_t>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("create_problem",
		&pda::implee2d::create_problem_for_pyB<ccumesh_t, euler_2d_t>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("create_problem",
		&pda::implee2d::create_problem_for_pyC<ccumesh_t, euler_2d_t>,
		pybind11::return_value_policy::take_ownership);

  // -----------------------
  // Swe 2d
  // -----------------------
  using swe_impl_2d_t = pda::implswe::EigenSwe2dApp<ccumesh_t>;

  using swe_2d_t = pressiodemoapps::PublicProblemMixinPy<swe_impl_2d_t>;
  pybind11::class_<swe_2d_t> swe2dProb(mTopLevel, "Swe2dProblem");
  pressiodemoappspy::impl::bindCommonApiMethods<swe_2d_t>(swe2dProb);
  swe2dProb.def("coriolis", &swe_2d_t::coriolis);
  swe2dProb.def("gravity", &swe_2d_t::gravity);

  mTopLevel.def("create_problem",
		&pda::implswe2d::create_problem_for_pyA<ccumesh_t, swe_2d_t>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("create_problem",
		&pda::implswe2d::create_problem_for_pyB<ccumesh_t, swe_2d_t>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("create_problem",
		&pda::implswe2d::create_problem_for_pyC<ccumesh_t, swe_2d_t>,
		pybind11::return_value_policy::take_ownership);

  // -----------------------
  // Euler 3d
  // -----------------------
  using euler_impl_3d_t = pda::ee::impl::EigenEuler3dApp<ccumesh_t>;

  using euler_3d_t = pressiodemoapps::PublicProblemMixinPy<euler_impl_3d_t>;
  pybind11::class_<euler_3d_t> ee3dProb(mTopLevel, "Euler3dProblem");
  pressiodemoappspy::impl::bindCommonApiMethods<euler_3d_t>(ee3dProb);
  ee3dProb.def("gamma", &euler_3d_t::gamma);

  mTopLevel.def("create_problem",
		&pda::implee3d::create_problem_for_pyA<ccumesh_t, euler_3d_t>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("create_problem",
		&pda::implee3d::create_problem_for_pyB<ccumesh_t, euler_3d_t>,
		pybind11::return_value_policy::take_ownership);

}
#endif

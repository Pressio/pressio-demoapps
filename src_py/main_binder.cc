
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
// #include "pressiodemoapps/swe2d.hpp"
// #include "pressiodemoapps/euler2d.hpp"
// #include "pressiodemoapps/euler3d.hpp"

namespace pressiodemoappspy{ namespace impl{

template<class cpp_t, class py_t>
void bindCommonExplicitApiMethods(py_t & appObj)
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

template<class cpp_t, class py_t>
void bindCommonImplicitApiMethods(py_t & appObj)
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

  appObj.def("createApplyJacobianResult", &cpp_t::createApplyJacobianResultRank1);
  appObj.def("createApplyJacobianResult", &cpp_t::createApplyJacobianResultRank2);
  appObj.def("applyJacobian", &cpp_t::applyJacobianRank1);
  appObj.def("applyJacobian", &cpp_t::applyJacobianRank2);
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

  pybind11::enum_<pda::Euler1d>(mParent, "Euler1d")
    .value("PeriodicSmooth", pda::Euler1d::PeriodicSmooth)
    .value("Sod",	     pda::Euler1d::Sod)
    .value("Lax",	     pda::Euler1d::Lax);

  // pybind11::enum_<pda::Euler2d>(mParent, "Euler2d")
  //   .value("PeriodicSmooth", pda::Euler2d::PeriodicSmooth)
  //   .value("KelvinHelmholtz", pda::Euler2d::KelvinHelmholtz)
  //   .value("SedovFull",	     pda::Euler2d::SedovFull)
  //   .value("SedovSymmetry",  pda::Euler2d::SedovSymmetry)
  //   .value("Riemann",	     pda::Euler2d::Riemann)
  //   .value("NormalShock",    pda::Euler2d::NormalShock)
  //   .value("DoubleMachReflection", pda::Euler2d::DoubleMachReflection);

  // pybind11::enum_<pda::Swe2d>(mParent, "Swe2d")
  //   .value("SlipWall", pda::Swe2d::SlipWall);

  // pybind11::enum_<pda::Euler3d>(mParent, "Euler3d")
  //   .value("PeriodicSmooth", pda::Euler3d::PeriodicSmooth)
  //   .value("SedovSymmetry",  pda::Euler3d::SedovSymmetry);
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
    mParent.def("loadCellCenterUniformMesh",
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
  using ad1d_t = pressiodemoapps::impladv::EigenAdvection1dApp<ccumesh_t>;

  using ad1d_exp_t = pressiodemoapps::ExplicitAdapterMixin<ad1d_t>;
  pybind11::class_<ad1d_exp_t> adv1dProbExp(mTopLevel, "Advection1dProblemExp");
  pressiodemoappspy::impl::bindCommonExplicitApiMethods<ad1d_exp_t>(adv1dProbExp);

  using ad1d_imp_t = pressiodemoapps::ImplicitAdapterMixinPy<ad1d_t>;
  pybind11::class_<ad1d_imp_t> adv1dProbImp(mTopLevel, "Advection1dProblemImp");
  pressiodemoappspy::impl::bindCommonImplicitApiMethods<ad1d_imp_t>(adv1dProbImp);

  mTopLevel.def("createExplicitProblem",
		&pda::impladv::createProblemForPy<ccumesh_t, ad1d_exp_t>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("createImplicitProblem",
		&pda::impladv::createProblemForPy<ccumesh_t, ad1d_imp_t>,
		pybind11::return_value_policy::take_ownership);

  // -----------------------
  // diffusion-reaction 1d
  // -----------------------
  using diffreac_1d_t = pda::impldiffreac::EigenDiffReac1dApp<ccumesh_t>;

  using diffreac_1d_exp_t = pressiodemoapps::ExplicitAdapterMixin<diffreac_1d_t>;
  pybind11::class_<diffreac_1d_exp_t> diffReac1dProbExp(mTopLevel, "DiffusionReaction1dProblemExp");
  pressiodemoappspy::impl::bindCommonExplicitApiMethods<diffreac_1d_exp_t>(diffReac1dProbExp);

  using diffreac_1d_imp_t = pressiodemoapps::ImplicitAdapterMixinPy<diffreac_1d_t>;
  pybind11::class_<diffreac_1d_imp_t> diffReac1dProbImp(mTopLevel, "DiffusionReaction1dProblemImp");
  pressiodemoappspy::impl::bindCommonImplicitApiMethods<diffreac_1d_imp_t>(diffReac1dProbImp);

  mTopLevel.def("createExplicitProblem",
		&pda::impldiffreac::createProblemForPyA<ccumesh_t, diffreac_1d_exp_t>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("createExplicitProblem",
		&pda::impldiffreac::createProblemForPyB<ccumesh_t, diffreac_1d_exp_t>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("createExplicitProblem",
		&pda::impldiffreac::createProblemForPyC<ccumesh_t, diffreac_1d_exp_t>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("createImplicitProblem",
		&pda::impldiffreac::createProblemForPyA<ccumesh_t, diffreac_1d_imp_t>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("createImplicitProblem",
		&pda::impldiffreac::createProblemForPyB<ccumesh_t, diffreac_1d_imp_t>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("createImplicitProblem",
		&pda::impldiffreac::createProblemForPyC<ccumesh_t, diffreac_1d_imp_t>,
		pybind11::return_value_policy::take_ownership);

  // -----------------------
  // Euler 1d
  // -----------------------
  using ee_1d_t = pda::implee1d::EigenEuler1dApp<ccumesh_t>;

  using ee_1d_exp_t = pressiodemoapps::ExplicitAdapterMixin<ee_1d_t>;
  pybind11::class_<ee_1d_exp_t> ee1dProbExp(mTopLevel, "Euler1dProblemExp");
  pressiodemoappspy::impl::bindCommonExplicitApiMethods<ee_1d_exp_t>(ee1dProbExp);

  using ee_1d_imp_t = pressiodemoapps::ImplicitAdapterMixinPy<ee_1d_t>;
  pybind11::class_<ee_1d_imp_t> ee1dProbImp(mTopLevel, "Euler1dProblemImp");
  pressiodemoappspy::impl::bindCommonImplicitApiMethods<ee_1d_imp_t>(ee1dProbImp);

  mTopLevel.def("createExplicitProblem",
		&pda::implee1d::createProblemForPyA<ccumesh_t, ee_1d_exp_t>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("createExplicitProblem",
		&pda::implee1d::createProblemForPyB<ccumesh_t, ee_1d_exp_t>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("createImplicitProblem",
		&pda::implee1d::createProblemForPyA<ccumesh_t, ee_1d_imp_t>,
		pybind11::return_value_policy::take_ownership);

  mTopLevel.def("createImplicitProblem",
		&pda::implee1d::createProblemForPyB<ccumesh_t, ee_1d_imp_t>,
		pybind11::return_value_policy::take_ownership);

  // // -----------------------
  // // Euler 2d
  // // -----------------------
  // using ee2d_t =
  //   pda::ee::impl::Euler2dAppRhsOnly<
  //     pressiodemoappspy::scalar_t,
  //   ccumesh_t,
  //   pressiodemoappspy::py_cstyle_arr_sc, // state type
  //   pressiodemoappspy::py_cstyle_arr_sc, // velo type
  //   pressiodemoappspy::py_cstyle_arr_sc  // ghost type
  //   >;

  // pybind11::class_<ee2d_t> ee2dClass(mTopLevel, "Euler2dProblem");
  // pressiodemoappspy::impl::bindCommonApiMethods<ee2d_t>(ee2dClass);
  // ee2dClass.def("gamma", &ee2d_t::gamma);

  // mTopLevel.def("createProblem",
  // 		&pda::createEuler2dForPyA<ccumesh_t, ee2d_t>,
  // 		pybind11::return_value_policy::take_ownership);

  // mTopLevel.def("createProblem",
  // 		&pda::createEuler2dForPyB<ccumesh_t, ee2d_t>,
  // 		pybind11::return_value_policy::take_ownership);

  // mTopLevel.def("createProblem",
  // 		&pda::createEuler2dForPyC<ccumesh_t, ee2d_t>,
  // 		pybind11::return_value_policy::take_ownership);


  // // -----------------------
  // // Swe 2d
  // // -----------------------
  // using swe2d_t =
  //   pda::implswe::Swe2dAppT<
  //     pressiodemoappspy::scalar_t,
  //   ccumesh_t,
  //   pressiodemoappspy::py_cstyle_arr_sc, // state type
  //   pressiodemoappspy::py_cstyle_arr_sc, // velo type
  //   pressiodemoappspy::py_cstyle_arr_sc  // ghost type
  //   >;

  // pybind11::class_<swe2d_t> swe2dClass(mTopLevel, "Swe2dProblem");
  // pressiodemoappspy::impl::bindCommonApiMethods<swe2d_t>(swe2dClass);
  // swe2dClass.def("coriolis", &swe2d_t::coriolis);
  // swe2dClass.def("gravity", &swe2d_t::gravity);

  // mTopLevel.def("createProblem",
  // 		&pda::createSwe2dForPyA<ccumesh_t, swe2d_t>,
  // 		pybind11::return_value_policy::take_ownership);

  // mTopLevel.def("createProblem",
  // 		&pda::createSwe2dForPyB<ccumesh_t, swe2d_t>,
  // 		pybind11::return_value_policy::take_ownership);

  // mTopLevel.def("createProblem",
  // 		&pda::createSwe2dForPyC<ccumesh_t, swe2d_t>,
  // 		pybind11::return_value_policy::take_ownership);

  // // -----------------------
  // // Euler 3d
  // // -----------------------
  // using ee3d_t =
  //   pda::ee::impl::Euler3dAppT<
  //     pressiodemoappspy::scalar_t,
  //   ccumesh_t,
  //   pressiodemoappspy::py_cstyle_arr_sc, // state type
  //   pressiodemoappspy::py_cstyle_arr_sc, // velo type
  //   pressiodemoappspy::py_cstyle_arr_sc  // ghost type
  //   >;

  // pybind11::class_<ee3d_t> ee3dClass(mTopLevel, "Euler3dProblem");
  // pressiodemoappspy::impl::bindCommonApiMethods<ee3d_t>(ee3dClass);
  // ee3dClass.def("gamma", &ee3d_t::gamma);

  // mTopLevel.def("createProblem",
  // 		&pda::createEuler3dForPyC<ccumesh_t, ee3d_t>,
  // 		pybind11::return_value_policy::take_ownership);
}
#endif

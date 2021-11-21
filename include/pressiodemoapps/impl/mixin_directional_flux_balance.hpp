
#ifndef PRESSIODEMOAPPS_MIXIN_DIRECTIONAL_FLUX_BALANCE_HPP_
#define PRESSIODEMOAPPS_MIXIN_DIRECTIONAL_FLUX_BALANCE_HPP_

namespace pressiodemoapps{ namespace impl{

template<class Parent, int ndpc, class ContainerType, class ScalarType>
struct ComputeDirectionalFluxBalance : Parent
{
private:
  ContainerType & m_V;
  ScalarType m_hInverse;

public:
  template<class ...Args>
  ComputeDirectionalFluxBalance(ContainerType & V,
				ScalarType hInverse,
				Args && ...args)
    : Parent(std::forward<Args>(args)...), m_V(V), m_hInverse(hInverse){}

  template<class index_t, class ...Args2>
  void operator()(index_t smPt, Args2 && ...args2)
  {
    Parent::operator()(smPt, std::forward<Args2>(args2)...);

    const auto & fluxL = Parent::fluxLeft();
    const auto & fluxR = Parent::fluxRight();
    const auto vIndexCurrentCellFirstDof = smPt*ndpc;
    for (int dof=0; dof<ndpc; ++dof){
      m_V(vIndexCurrentCellFirstDof+dof) += m_hInverse*(fluxL(dof) - fluxR(dof));
    }
  }
};

template<class Parent, class ContainerType, class ScalarType>
struct ComputeDirectionalFluxBalance<Parent, 1, ContainerType, ScalarType> : Parent
{
private:
  ContainerType & m_V;
  ScalarType m_hInverse;

public:
  template<class ...Args>
  ComputeDirectionalFluxBalance(ContainerType & V,
				ScalarType hInverse,
				Args && ...args)
    : Parent(std::forward<Args>(args)...), m_V(V), m_hInverse(hInverse){}

  template<class index_t, class ...Args2>
  void operator()(index_t smPt, Args2 && ...args2)
  {
    Parent::operator()(smPt, std::forward<Args2>(args2)...);

    const auto & fluxL = Parent::fluxLeft();
    const auto & fluxR = Parent::fluxRight();
    m_V(smPt) += m_hInverse*(fluxL - fluxR);
  }
};

}} //end namespaces
#endif

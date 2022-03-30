
#ifndef PRESSIODEMOAPPS_MIXIN_DIRECTIONAL_FLUX_BALANCE_NON_TEMP_HPP_
#define PRESSIODEMOAPPS_MIXIN_DIRECTIONAL_FLUX_BALANCE_NON_TEMP_HPP_

namespace pressiodemoapps{ namespace impl{

template<class Parent, class ContainerType, class ScalarType>
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
  void operator()(index_t smPt, int ndpc, Args2 && ...args2)
  {
    Parent::operator()(smPt, ndpc, std::forward<Args2>(args2)...);

    const auto & fluxL = Parent::fluxLeft();
    const auto & fluxR = Parent::fluxRight();
    if (ndpc == 1){
      m_V(smPt) += m_hInverse*(fluxL(0) - fluxR(0));
    }
    else{
      const auto vIndexCurrentCellFirstDof = smPt*ndpc;
      for (int dof=0; dof<ndpc; ++dof){
	m_V(vIndexCurrentCellFirstDof+dof) += m_hInverse*(fluxL(dof) - fluxR(dof));
      }
    }
  }
};

}} //end namespaces
#endif

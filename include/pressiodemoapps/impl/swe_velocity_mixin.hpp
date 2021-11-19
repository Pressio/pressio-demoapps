
#ifndef PRESSIODEMOAPPS_SWE_MIXIN_CELL_VELOCITY_HPP_
#define PRESSIODEMOAPPS_SWE_MIXIN_CELL_VELOCITY_HPP_

namespace pressiodemoapps{ namespace implswe{

template<class Parent, class VeloType, class ScalarType>
struct CellVelocityNoForcing : Parent
{
private:
  VeloType & m_V;
  ScalarType m_hInv;

public:
  template<class ...Args>
  CellVelocityNoForcing(VeloType & V,
	       ScalarType hInv,
	       Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_V(V), m_hInv(hInv){}

  template<class index_t, class ...Args2>
  void operator()(index_t smPt, Args2 && ...args2)
  {
    Parent::operator()(smPt, std::forward<Args2>(args2)...);

    const auto & fluxL = Parent::fluxLeft();
    const auto & fluxR = Parent::fluxRight();
    const auto vIndex = smPt*3;
    m_V(vIndex)   += m_hInv*(fluxL(0) - fluxR(0));
    m_V(vIndex+1) += m_hInv*(fluxL(1) - fluxR(1));
    m_V(vIndex+2) += m_hInv*(fluxL(2) - fluxR(2));
  }
};

}} //end namespaces
#endif

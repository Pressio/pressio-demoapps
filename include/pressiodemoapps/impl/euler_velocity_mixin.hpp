
#ifndef PRESSIODEMOAPPS_EULER_MIXIN_CELL_VELOCITY_HPP_
#define PRESSIODEMOAPPS_EULER_MIXIN_CELL_VELOCITY_HPP_

namespace pressiodemoapps{ namespace ee{ namespace impl{

template<class Parent, int ndpc, class VeloType, class ScalarType>
struct CellVelocity : Parent
{
private:
  VeloType & m_V;
  ScalarType m_hInv;

public:
  template<class ...Args>
  CellVelocity(VeloType & V, ScalarType hInv, Args && ...args)
    : Parent(std::forward<Args>(args)...), m_V(V), m_hInv(hInv){}

  template<class index_t, class ...Args2>
  void operator()(index_t smPt, Args2 && ...args2)
  {
    Parent::operator()(smPt, std::forward<Args2>(args2)...);

    const auto & fluxL = Parent::fluxLeft();
    const auto & fluxR = Parent::fluxRight();
    const auto vIndexCurrentCellFirstDof = smPt*ndpc;
    for (int dof=0; dof<ndpc; ++dof){
      m_V(vIndexCurrentCellFirstDof+dof) += m_hInv*(fluxL(dof) - fluxR(dof));
    }
  }
};

}}} //end namespaces
#endif

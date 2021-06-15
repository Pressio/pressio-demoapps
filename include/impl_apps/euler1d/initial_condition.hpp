
#ifndef PRESSIODEMOAPPS_EE1D_IC_HPP_
#define PRESSIODEMOAPPS_EE1D_IC_HPP_

namespace pressiodemoapps{ namespace ee{

template<class state_type, class mesh_t, class scalar_type>
void euler1dsineInitialCondition(state_type & state,
				 const mesh_t & meshObj,
				 const scalar_type gamma)
{

  constexpr int numDofPerCell = 3;
  const auto & x= meshObj.viewX();
  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      prim[0] = 1. + 0.2*std::sin(M_PI*x(i));
      prim[1] = 1.;
      prim[2] = 1.;

      const auto ind = i*numDofPerCell;
      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = computeEnergyFromPrimitive(gamma, prim);
    }
}

template<class state_type, class mesh_t, class scalar_type>
void sod1dInitialCondition(state_type & state,
			   const mesh_t & meshObj,
			   const scalar_type gamma)
{

  constexpr int numDofPerCell = 3;
  const auto & x= meshObj.viewX();
  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      if (x(i) <= 0.0){
	prim[0] = 1.0;
	prim[1] = 0.0;
	prim[2] = 1.0;
      }

      if (x(i) > 0.0){
	prim[0] = 0.125;
	prim[1] = 0.0;
	prim[2] = 0.1;
      }

      const auto ind = i*numDofPerCell;
      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = computeEnergyFromPrimitive(gamma, prim);
    }
}

template<class state_type, class mesh_t, class scalar_type>
void lax1dInitialCondition(state_type & state,
			   const mesh_t & meshObj,
			   const scalar_type gamma)
{

  constexpr int numDofPerCell = 3;
  const auto & x= meshObj.viewX();
  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      const auto ind = i*numDofPerCell;
      if (x(i) <= 0.0){
	prim[0] = 0.445;
	prim[1] = 0.698;
	prim[2] = 3.528;
      }
      else if (x(i) > 0.0){
	prim[0] = 0.5;
	prim[1] = 0.;
	prim[2] = 0.571;
      }

      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = computeEnergyFromPrimitive(gamma, prim);
    }
}

}}//end namespace

#endif

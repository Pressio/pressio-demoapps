
#ifndef PRESSIODEMOAPPS_EE1D_IC_HPP_
#define PRESSIODEMOAPPS_EE1D_IC_HPP_

namespace pressiodemoapps{ namespace impl{

template<class state_type, class mesh_t, class scalar_type>
void euler1dsineInitialCondition(state_type & state,
				 const mesh_t & meshObj,
				 const scalar_type gamma)
{
  constexpr int numDofPerCell = 3;
  constexpr auto one  = static_cast<scalar_type>(1);

  const auto & x= meshObj.viewX();
  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      prim[0] = one + static_cast<scalar_type>(0.2)*std::sin(M_PI*x(i));
      prim[1] = one;
      prim[2] = one;

      const auto ind = i*numDofPerCell;
      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = eulerEquationsComputeEnergyFromPrimitive(gamma, prim);
    }
}

template<class state_type, class mesh_t, class scalar_type>
void sod1dInitialCondition(state_type & state,
			   const mesh_t & meshObj,
			   const scalar_type gamma)
{

  constexpr int numDofPerCell = 3;
  constexpr auto zero = static_cast<scalar_type>(0);
  constexpr auto one  = static_cast<scalar_type>(1);

  const auto & x= meshObj.viewX();
  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      if (x(i) <= zero){
	prim[0] = one;
	prim[1] = zero;
	prim[2] = one;
      }

      if (x(i) > zero){
	prim[0] = static_cast<scalar_type>(0.125);
	prim[1] = zero;
	prim[2] = static_cast<scalar_type>(0.1);
      }

      const auto ind = i*numDofPerCell;
      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = eulerEquationsComputeEnergyFromPrimitive(gamma, prim);
    }
}

template<class state_type, class mesh_t, class scalar_type>
void lax1dInitialCondition(state_type & state,
			   const mesh_t & meshObj,
			   const scalar_type gamma)
{

  constexpr int numDofPerCell = 3;
  constexpr auto zero = static_cast<scalar_type>(0);
  constexpr auto one  = static_cast<scalar_type>(1);

  const auto & x= meshObj.viewX();
  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      const auto ind = i*numDofPerCell;
      if (x(i) <= zero){
	prim[0] = static_cast<scalar_type>(0.445);
	prim[1] = static_cast<scalar_type>(0.698);
	prim[2] = static_cast<scalar_type>(3.528);
      }
      else if (x(i) > zero){
	prim[0] = static_cast<scalar_type>(0.5);
	prim[1] = zero;
	prim[2] = static_cast<scalar_type>(0.571);
      }

      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = eulerEquationsComputeEnergyFromPrimitive(gamma, prim);
    }
}

template<class state_type, class mesh_t, class scalar_type>
void shuOsherInitialCondition(state_type & state,
			   const mesh_t & meshObj,
			   const scalar_type gamma)
{

  constexpr int numDofPerCell = 3;
  constexpr auto zero    = static_cast<scalar_type>(0);
  constexpr auto one     = static_cast<scalar_type>(1);
  constexpr auto three   = static_cast<scalar_type>(3);
  constexpr auto four    = static_cast<scalar_type>(4);
  constexpr auto negFour = -four;
  constexpr auto five    = static_cast<scalar_type>(5);
  constexpr auto seven   = static_cast<scalar_type>(7);
  constexpr auto twentySeven  = static_cast<scalar_type>(27);
  constexpr auto thirtyOne  = static_cast<scalar_type>(31);

  const auto & x= meshObj.viewX();
  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      const auto ind = i*numDofPerCell;
      if (x(i) <= negFour){
	prim[0] = twentySeven/seven;
	prim[1] = static_cast<scalar_type>(2.629369);
	prim[2] = thirtyOne/three;
      }
      else{
	prim[0] = one + (one/five)*std::sin(five*x(i));
	prim[1] = zero;
	prim[2] = one;
      }

      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = eulerEquationsComputeEnergyFromPrimitive(gamma, prim);
    }
}

}}//end namespace
#endif

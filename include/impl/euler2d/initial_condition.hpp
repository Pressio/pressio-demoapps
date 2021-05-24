
#ifndef PRESSIODEMOAPPS_IC_HPP_
#define PRESSIODEMOAPPS_IC_HPP_

namespace pressiodemoapps{ namespace ee{

/*
  probid = 1: Sedov,
https://www.researchgate.net/publication/260967068_GENASIS_General_Astrophysical_Simulation_System_I_Refinable_Mesh_and_Nonrelativistic_Hydrodynamics
*/

template<class state_type, class mesh_t, class scalar_type>
void sedov2dInitialCondition(state_type & state,
			     const mesh_t & meshObj,
			     const int numDofPerCell,
			     const scalar_type gamma)
{

  const auto dx  = meshObj.dx();
  const auto dy  = meshObj.dy();

  const auto sRad = 2.*std::min(dx, dy);
  const auto x= meshObj.viewX();
  const auto y= meshObj.viewY();
  const auto gammaMinusOne = gamma - 1.;
  const auto gammaMinusOneInv = 1./gammaMinusOne;

  std::array<scalar_type, 4> prim;
  for (int i=0; i<x.size(); ++i)
    {

      const auto ind = i*numDofPerCell;
      const auto diffx = x(i);
      const auto diffy = y(i);
      const auto myR = std::sqrt(diffx*diffx + diffy*diffy);

      if (myR <= sRad)
	{
	  prim[0] = 1.0;
	  prim[1] = 0.0;
	  prim[2] = 0.0;
	  prim[3] = (gammaMinusOne)/(M_PI*sRad*sRad);

	  state(ind)   = prim[0];
	  state(ind+1) = prim[0]*prim[1];
	  state(ind+2) = prim[0]*prim[2];
	  state(ind+3) = computeEnergyFromPrimitive2(gammaMinusOneInv, prim);
	}
      else
	{
	  prim[0] = 1.;
	  prim[1] = 0.0;
	  prim[2] = 0.0;
	  prim[3] = 5.e-5;

	  state(ind)   = prim[0];
	  state(ind+1) = prim[0]*prim[1];
	  state(ind+2) = prim[0]*prim[2];
	  state(ind+3) = computeEnergyFromPrimitive2(gammaMinusOneInv, prim);
	}
    }

}


/*
  first case page 15 of
  https://www.researchgate.net/publication/269636534_A_Compact_Third-Order_Gas-Kinetic_Scheme_for_Compressible_Euler_and_Navier-Stokes_Equations
 */

template<class state_type, class mesh_t, class scalar_type>
void riemann2dInitialCondition1(state_type & state,
				const mesh_t & meshObj,
				const int numDofPerCell,
				const scalar_type gamma)
{

  const auto dx  = meshObj.dx();
  const auto dy  = meshObj.dy();
  const auto x= meshObj.viewX();
  const auto y= meshObj.viewY();
  const auto gammaMinusOne = gamma - 1.;
  const auto gammaMinusOneInv = 1./gammaMinusOne;

  constexpr auto x0 = 0.5;
  constexpr auto y0 = 0.5;

  std::array<scalar_type, 4> prim;
  for (int i=0; i<x.size(); ++i)
    {

      if (x[i] >= x0 and y[i] >= y0){
	prim = {0.5313, 0., 0., 0.4};
      }
      else if (x[i] < x0 and y[i] >= y0){
	prim = {1., 0.7276, 0., 1.};
      }
      else if (x[i] < x0 and y[i] < y0){
	prim = {0.8, 0., 0., 1.};
      }
      else if (x[i] > x0 and y[i] < y0){
	prim = {1., 0., 0.7276, 1.};
      }

      const auto ind = i*numDofPerCell;
      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = prim[0]*prim[2];
      state(ind+3) = computeEnergyFromPrimitive2(gammaMinusOneInv, prim);
   }

}


/*
  IC2:
  http://www.amsc-ouc.ac.cn/Files/Papers/2016_Don_Hybrid%20Compact-WENO%20finite%20difference%20scheme%20with%20conjugate%20Fourier%20shock%20detection%20algorithm%20for%20hyperbolic%20conservation%20laws.pdf
 */
template<class state_type, class mesh_t, class scalar_type>
void riemann2dInitialCondition2(state_type & state,
				const mesh_t & meshObj,
				const int numDofPerCell,
				const scalar_type gamma)
{

  const auto dx  = meshObj.dx();
  const auto dy  = meshObj.dy();
  const auto x0 = 0.8;
  const auto y0 = 0.8;
  const auto gammaMinusOne = gamma - 1.;
  const auto gammaMinusOneInv = 1./gammaMinusOne;

  const auto x= meshObj.viewX();
  const auto y= meshObj.viewY();
  std::array<scalar_type, 4> prim;
  for (int i=0; i<x.size(); ++i)
    {
      if (x[i] >= x0 and y[i] >= y0){
	prim = {1.5, 0., 0., 1.5};
      }
      else if (x[i] < x0 and y[i] >= y0){
	prim = {0.5323, 1.206, 0., 0.3};
      }
      else if (x[i] < x0 and y[i] < y0){
	prim = {0.138, 1.206, 1.206, 0.029};
      }
      else if (x[i] > x0 and y[i] < y0){
	prim = {0.5323, 0., 1.206, 0.3};
      }

      const auto ind = i*numDofPerCell;
      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = prim[0]*prim[2];
      state(ind+3) = computeEnergyFromPrimitive2(gammaMinusOneInv, prim);
    }
}

}}//end namespace

#endif

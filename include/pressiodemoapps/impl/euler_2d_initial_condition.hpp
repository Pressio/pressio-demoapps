
#ifndef PRESSIODEMOAPPS_EULER2D_IC_HPP_
#define PRESSIODEMOAPPS_EULER2D_IC_HPP_

namespace pressiodemoapps{
namespace impleuler2d{

template<class state_type, class mesh_t, class scalar_type>
void sin2dEulerIC(state_type & state,
		  const mesh_t & meshObj,
		  const scalar_type gamma)
{
  constexpr int numDofPerCell = 4;
  constexpr auto one  = static_cast<scalar_type>(1);
  constexpr auto five = static_cast<scalar_type>(5);

  const auto &x= meshObj.viewX();
  const auto &y= meshObj.viewY();
  const auto gammaMinusOne = gamma - one;
  const auto gammaMinusOneInv = one/gammaMinusOne;

  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      const auto ind = i*numDofPerCell;

      prim[0] = one + (one/five)*std::sin(M_PI*(x(i)+y(i)));
      prim[1] = one;
      prim[2] = one;
      prim[3] = one;
      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = prim[0]*prim[2];
      state(ind+3) = eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, prim);
    }
}

template<class state_type, class mesh_t, class scalar_type>
void KelvinHelmholtzIC(state_type & state,
      const mesh_t & meshObj,
      const scalar_type gamma)
{
  constexpr int numDofPerCell = 4;
  constexpr auto one  = static_cast<scalar_type>(1);

  const auto &x= meshObj.viewX();
  const auto &y= meshObj.viewY();
  const auto gammaMinusOne = gamma - one;
  const auto gammaMinusOneInv = one/gammaMinusOne;

  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      const auto ind = i*numDofPerCell;
      scalar_type freq = 4.;
      scalar_type mag = 0.025;
      scalar_type pert = mag*std::cos(2.*3.14159265/10.*freq*x(i));
      if (y(i) > -2 + pert && y(i) < 2 + pert){
        state(ind) = 2.;
        state(ind+1) = state(ind)*0.5;
        state(ind+2) = 0.;
        state(ind+3) = 2.5 / (gammaMinusOne) + 0.5/state(ind)*(state(ind+1)*state(ind+1) + state(ind+2)*state(ind+2) );
        //ainf = sqrt(gamma p / rho) = 1.3228756555322954
        //Minf = 0.3779644730092272
      }
      else{
        state(ind) = 1.;
        state(ind+1) = -state(ind)*0.5;
        state(ind+2) = 0.;
        state(ind+3) = 2.5 / (gammaMinusOne) + 0.5/state(ind)*(state(ind+1)*state(ind+1) + state(ind+2)*state(ind+2) );
        //ainf = np.sqrt(gamma * p / rho) = 1.8708286933869707
        //Minf = 0.2672612419124244
      }
    }
}


/*
  probid = 1: Sedov,
  https://www.researchgate.net/publication/260967068_GENASIS_General_Astrophysical_Simulation_System_I_Refinable_Mesh_and_Nonrelativistic_Hydrodynamics
*/

template<class state_type, class mesh_t, class scalar_type>
void sedov2dIC(state_type & state,
	       const mesh_t & meshObj,
	       const scalar_type gamma)
{
  constexpr int numDofPerCell = 4;
  constexpr auto zero = static_cast<scalar_type>(0);
  constexpr auto one  = static_cast<scalar_type>(1);
  constexpr auto two  = static_cast<scalar_type>(2);

  const auto dx  = meshObj.dx();
  const auto dy  = meshObj.dy();
  const auto &x= meshObj.viewX();
  const auto &y= meshObj.viewY();
  const auto gammaMinusOne = gamma - one;
  const auto gammaMinusOneInv = one/gammaMinusOne;
  const auto sRad = two*std::min(dx, dy);

  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {

      const auto ind = i*numDofPerCell;
      const auto diffx = x(i);
      const auto diffy = y(i);
      const auto myR = std::sqrt(diffx*diffx + diffy*diffy);

      if (myR <= sRad)
	{
	  prim[0] = one;
	  prim[1] = zero;
	  prim[2] = zero;
	  prim[3] = (gammaMinusOne)/(M_PI*sRad*sRad);

	  state(ind)   = prim[0];
	  state(ind+1) = prim[0]*prim[1];
	  state(ind+2) = prim[0]*prim[2];
	  state(ind+3) = eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, prim);
	}

      else
	{
	  prim[0] = one;
	  prim[1] = zero;
	  prim[2] = zero;
	  prim[3] = static_cast<scalar_type>(5.e-5);

	  state(ind)   = prim[0];
	  state(ind+1) = prim[0]*prim[1];
	  state(ind+2) = prim[0]*prim[2];
	  state(ind+3) = eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, prim);
	}
    }
}

/*
the 3 for the radius is taken from
http://flash.uchicago.edu/site/flashcode/user_support/flash_ug_devel/node184.html#SECTION010114000000000000000

energy value taken from
https://reader.elsevier.com/reader/sd/pii/S002199911400477X?token=658F08D28B5C7A6A97E6F4478FD494699F3C8DF23970A256F06E501B7B136F9A6A540EEA749F28AC2AF4A6A7993A8517&originRegion=eu-west-1&originCreation=20210611123033

this problem is very delicate, so the IC has been tried many times.
Seems to be sort of magic seting up these numbers...
*/

template<class state_type, class mesh_t, class scalar_type>
void sedov2dsymmetryIC(state_type & state,
		       const mesh_t & meshObj,
		       const scalar_type gamma)
{
  constexpr int numDofPerCell = 4;
  constexpr auto zero  = static_cast<scalar_type>(0);
  constexpr auto one   = static_cast<scalar_type>(1);

  const auto dx  = meshObj.dx();
  const auto dy  = meshObj.dy();
  const auto &x= meshObj.viewX();
  const auto &y= meshObj.viewY();
  const auto gammaMinusOne = gamma - one;
  const auto gammaMinusOneInv = one/gammaMinusOne;
  const auto sRad = 3.*std::min(dx, dy);

  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {

      const auto ind = i*numDofPerCell;
      const auto diffx = x(i);
      const auto diffy = y(i);
      const auto myR = std::sqrt(diffx*diffx + diffy*diffy);

      if (myR <= sRad)
	{
	  prim[0] = one;
	  prim[1] = zero;
	  prim[2] = zero;
	  prim[3] = gammaMinusOne*0.851072/(M_PI*sRad*sRad);

	  state(ind)   = prim[0];
	  state(ind+1) = prim[0]*prim[1];
	  state(ind+2) = prim[0]*prim[2];
	  state(ind+3) = eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, prim);
	}

      else
	{
	  prim[0] = one;
	  prim[1] = zero;
	  prim[2] = zero;
	  prim[3] = static_cast<scalar_type>(2.5e-5);

	  state(ind)   = prim[0];
	  state(ind+1) = prim[0]*prim[1];
	  state(ind+2) = prim[0]*prim[2];
	  state(ind+3) = eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, prim);
	}
    }
}

/*
  first case page 15 of
  https://www.researchgate.net/publication/269636534_A_Compact_Third-Order_Gas-Kinetic_Scheme_for_Compressible_Euler_and_Navier-Stokes_Equations
*/

template<class state_type, class mesh_t, class scalar_type>
void riemann2dIC1(state_type & state,
		  const mesh_t & meshObj,
		  const scalar_type gamma)
{
  constexpr int numDofPerCell = 4;
  constexpr auto zero = static_cast<scalar_type>(0);
  constexpr auto one = static_cast<scalar_type>(1);
  constexpr auto two = static_cast<scalar_type>(2);

  const auto & x= meshObj.viewX();
  const auto & y= meshObj.viewY();
  const auto gammaMinusOne = gamma - one;
  const auto gammaMinusOneInv = one/gammaMinusOne;
  constexpr auto x0 = one/two;
  constexpr auto y0 = one/two;

  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {

      if (x(i) >= x0 and y(i) >= y0){
	prim = {0.5313, zero, zero, 0.4};
      }
      else if (x(i) < x0 and y(i) >= y0){
	prim = {one, 0.7276, zero, one};
      }
      else if (x(i) < x0 and y(i) < y0){
	prim = {0.8, zero, zero, one};
      }
      else if (x(i) > x0 and y(i) < y0){
	prim = {one, zero, 0.7276, one};
      }

      const auto ind = i*numDofPerCell;
      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = prim[0]*prim[2];
      state(ind+3) = eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, prim);
    }

}


/*
  IC2:
  http://www.amsc-ouc.ac.cn/Files/Papers/2016_Don_Hybrid%20Compact-WENO%20finite%20difference%20scheme%20with%20conjugate%20Fourier%20shock%20detection%20algorithm%20for%20hyperbolic%20conservation%20laws.pdf
*/
template<class state_type, class mesh_t, class scalar_type>
void riemann2dIC2(state_type & state,
		  const mesh_t & meshObj,
		  const scalar_type gamma)
{
  constexpr int numDofPerCell = 4;
  constexpr auto zero = static_cast<scalar_type>(0);
  constexpr auto one = static_cast<scalar_type>(1);

  const auto x0 = static_cast<scalar_type>(0.8);
  const auto y0 = static_cast<scalar_type>(0.8);
  const auto gammaMinusOne = gamma - one;
  const auto gammaMinusOneInv = one/gammaMinusOne;

  const auto & x= meshObj.viewX();
  const auto & y= meshObj.viewY();
  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      if (x(i) >= x0 and y(i) >= y0){
	prim = {1.5, zero, zero, 1.5};
      }
      else if (x(i) < x0 and y(i) >= y0){
	prim = {0.5323, 1.206, zero, 0.3};
      }
      else if (x(i) < x0 and y(i) < y0){
	prim = {0.138, 1.206, 1.206, 0.029};
      }
      else if (x(i) > x0 and y(i) < y0){
	prim = {0.5323, zero, 1.206, 0.3};
      }

      const auto ind = i*numDofPerCell;
      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = prim[0]*prim[2];
      state(ind+3) = eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, prim);
    }
}

/*
  normal shock propagating to right
*/
template<class state_type, class mesh_t, class scalar_type>
void normalShock2dIC(state_type & state,
		     const mesh_t & meshObj,
		     const scalar_type gamma)
{
  constexpr int numDofPerCell = 4;
  constexpr auto zero = static_cast<scalar_type>(0);
  constexpr auto one  = static_cast<scalar_type>(1);
  constexpr auto six  = static_cast<scalar_type>(6);

  constexpr scalar_type machShock{9};
  constexpr auto angle = zero;
  constexpr auto xShock = one/six;

  const auto & x  = meshObj.viewX();
  const auto gammaMinusOne = gamma - one;
  const auto gammaMinusOneInv = one/gammaMinusOne;

  // Compute Shock conditions
  std::array<scalar_type, numDofPerCell> primPreShock;
  std::array<scalar_type, numDofPerCell> primPostShock;

  // Density, shock normal velocity, pressure
  primPreShock[0] = gamma;
  primPreShock[1] = zero;
  primPreShock[2] = zero;
  primPreShock[3] = one;
  ::pressiodemoapps::ee::computePostShockConditionsFromPreshockAtRest(primPostShock,
								      primPreShock,
								      angle,
								      machShock,
								      gamma);

  for (int i=0; i<x.size(); ++i)
    {
      const auto ind = i*numDofPerCell;
      if (x(i) < xShock){
	state(ind)   = primPostShock[0];
	state(ind+1) = primPostShock[0]*primPostShock[1];
	state(ind+2) = primPostShock[0]*primPostShock[2];
	state(ind+3) = eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, primPostShock);
      }
      else {
        state(ind)   = primPreShock[0];
        state(ind+1) = primPreShock[0]*primPreShock[1];
        state(ind+2) = primPreShock[0]*primPreShock[2];
        state(ind+3) = eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, primPreShock);
      }
    }
}


/*
  Double Mach Reflection Wave,
  https://homepages.dias.ie/jmackey/jmac/node10.html
*/
template<class state_type, class mesh_t, class scalar_type>
void doubleMachReflection2dIC(state_type & state,
			      const mesh_t & meshObj,
			      const scalar_type gamma)
{
  constexpr int numDofPerCell = 4;
  constexpr auto zero = static_cast<scalar_type>(0);
  constexpr auto one  = static_cast<scalar_type>(1);
  constexpr auto six  = static_cast<scalar_type>(6);

  constexpr auto xShockBottom = one/six;
  constexpr scalar_type machShock{10};
  constexpr scalar_type angle{M_PI/six};
  const auto shockSlope = std::tan(angle);

  const auto & x  = meshObj.viewX();
  const auto & y  = meshObj.viewY();
  const auto gammaMinusOne = gamma - one;
  const auto gammaMinusOneInv = one/gammaMinusOne;

  // Compute Shock conditions
  std::array<scalar_type, numDofPerCell> primPreShock;
  std::array<scalar_type, numDofPerCell> primPostShock;

  // Density, shock normal velocity, pressure
  primPreShock[0] = gamma;
  primPreShock[1] = zero;
  primPreShock[2] = zero;
  primPreShock[3] = one;
  ::pressiodemoapps::ee::computePostShockConditionsFromPreshockAtRest(primPostShock,
								      primPreShock,
								      -angle,
								      machShock,
								      gamma);

  for (int i=0; i<x.size(); ++i)
    {
      const auto ind = i*numDofPerCell;
      const auto xShock = xShockBottom + shockSlope * y(i);

      if (x(i) < xShock)
	{
	  state(ind)   = primPostShock[0];
	  state(ind+1) = primPostShock[0]*primPostShock[1];
	  state(ind+2) = primPostShock[0]*primPostShock[2];
	  state(ind+3) = eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, primPostShock);
	}

      else {
        state(ind)   = primPreShock[0];
        state(ind+1) = primPreShock[0]*primPreShock[1];
        state(ind+2) = primPreShock[0]*primPreShock[2];
        state(ind+3) = eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, primPreShock);
      }

    }
}

}}//end namespace
#endif

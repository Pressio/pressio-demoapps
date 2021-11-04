
#ifndef PRESSIODEMOAPPS_EULER3D_IC_HPP_
#define PRESSIODEMOAPPS_EULER3D_IC_HPP_

namespace pressiodemoapps{ namespace ee{

template<class state_type, class mesh_t, class scalar_type>
void euler3dsmoothInitialCondition(state_type & state,
				   const mesh_t & meshObj,
				   const scalar_type gamma)
{
  constexpr int numDofPerCell = 5;
  const auto &x= meshObj.viewX();
  const auto &y= meshObj.viewY();
  const auto &z= meshObj.viewZ();
  const auto gammaMinusOne = gamma - 1.;
  const auto gammaMinusOneInv = 1./gammaMinusOne;

  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {

      const auto ind = i*numDofPerCell;
      prim[0] = 1.0 + 0.2*std::sin(M_PI*(x(i)+y(i)+z(i)));
      prim[1] = 1.0;
      prim[2] = 1.0;
      prim[3] = 1.0;
      prim[4] = 1.;
      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = prim[0]*prim[2];
      state(ind+3) = prim[0]*prim[3];
      state(ind+4) = eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, prim);
    }
}

template<class state_type, class mesh_t, class scalar_type>
void sedov3dInitialCondition(state_type & state,
			     const mesh_t & meshObj,
			     const scalar_type gamma)
{
  constexpr int numDofPerCell = 5;

  const auto dx  = meshObj.dx();
  const auto dy  = meshObj.dy();
  const auto dz  = meshObj.dz();
  const auto sRad = 3.*std::min(dx, std::min(dy, dz) );

  const auto &x= meshObj.viewX();
  const auto &y= meshObj.viewY();
  const auto &z= meshObj.viewZ();
  const auto gammaMinusOne = gamma - 1.;
  const auto gammaMinusOneInv = 1./gammaMinusOne;

  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      const auto ind = i*numDofPerCell;
      const auto distX = x(i);
      const auto distY = y(i);
      const auto distZ = z(i);
      const auto xsq = distX*distX;
      const auto ysq = distY*distY;
      const auto zsq = distZ*distZ;
      const auto myR = std::sqrt(xsq+ysq+zsq);

      if (myR <= sRad)
	{
	  prim[0] = 1.0;
	  prim[1] = 0.0;
	  prim[2] = 0.0;
	  prim[3] = 0.0;
	  prim[4] = (3.*gammaMinusOne*0.851072)/(4.*M_PI*sRad*sRad*sRad);
	}

      else
	{
	  prim[0] = 1.;
	  prim[1] = 0.0;
	  prim[2] = 0.0;
	  prim[3] = 0.0;
	  prim[4] = 2.5e-5;
	}

      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = prim[0]*prim[2];
      state(ind+3) = prim[0]*prim[3];
      state(ind+4) = eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, prim);
    }
}

}}//end namespace

#endif

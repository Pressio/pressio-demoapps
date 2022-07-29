
#include "pressiodemoapps/weno.hpp"
#include <array>
#include <iostream>
#include <iomanip>
#include <vector>

bool api1()
{
  using scalar_t = double;
  using vec_t = std::array<scalar_t, 7>;
  vec_t q{8., 5., 0., 1., 2., 1., 4.};

  std::array<scalar_t, 4> f = {};
  pressiodemoapps::weno5(f[0], f[1], f[2], f[3], q);

  for (auto & it : f){
    std::cout << std::setprecision(15) << it << " \n";
  }

  std::array<scalar_t, 4> goldF =
    {0.498169371091894,
     0.498239224509689,
     1.50249211669906,
     1.53308488724067};

  for (int i=0; i<4; ++i){
    const auto err = std::abs(goldF[i] - f[i]);
    if (err > 1e-13){
      return false;
    }
  }
  return true;
}

bool api2()
{
  using scalar_t = double;
  using vec_t = std::vector<scalar_t>;
  vec_t q({8.,8., 5.,5., 0.,0., 1.,1., 2.,2., 1.,1., 4.,4.});

  vec_t fLn = {0.,0.};
  vec_t fLp = {0.,0.};
  vec_t fRn = {0.,0.};
  vec_t fRp = {0.,0.};

  pressiodemoapps::weno5(fLn, fLp, fRn, fRp, q,
			 0, 2, 4, 8, 10, 12, 6, 2);

  std::array<scalar_t, 4> goldF =
    {0.498169371091894,
     0.498239224509689,
     1.50249211669906,
     1.53308488724067};

  for (int k=0; k<2; ++k)
  {
    const auto e1 = std::abs(goldF[0] - fLn[k]);
    if (e1 > 1e-13){ return false; }

    const auto e2 = std::abs(goldF[1] - fLp[k]);
    if (e2 > 1e-13){ return false; }

    const auto e3 = std::abs(goldF[2] - fRn[k]);
    if (e3 > 1e-13){ return false; }

    const auto e4 = std::abs(goldF[3] - fRp[k]);
    if (e4 > 1e-13){ return false; }

  }
  return true;
}

int main()
{
  const auto s1 = api1();
  const auto s2 = api2();
  if(s1 and s2){
    std::puts("PASS");
  }
  else{
    std::puts("FAILED");
  }

  return 0;
}

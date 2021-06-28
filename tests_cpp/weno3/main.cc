
#include "weno.hpp"
#include <array>
#include <iostream>
#include <iomanip>
#include <vector>

using scalar_t = double;

bool api1()
{
  std::array<scalar_t, 4> goldF =
  { 1.85840705145274,
    -0.451807145304084,
    0.487730039670285,
    0.5
  };

  using vec_t = std::array<scalar_t, 5>;
  vec_t q{5., 3., 0., 1., 2.};

  std::array<scalar_t, 4> f = {};
  pressiodemoapps::weno3(f[0], f[1], f[2], f[3], q);

  for (auto & it : f){
    std::cout << std::setprecision(15) << it << " \n";
  }

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
  std::array<scalar_t, 4> goldF0 =
  { 1.85840705145274,
    -0.451807145304084,
    0.487730039670285,
    0.5
  };

  std::array<scalar_t, 4> goldF1 =
  {
    2.85840705145274,
    0.548192854695916,
    1.48773003967029,
    1.5,
  };

  using vec_t = std::vector<scalar_t>;
  vec_t q({5.,6., 3.,4., 0.,1., 1.,2., 2.,3.});

  vec_t fLn = {0.,0.};
  vec_t fLp = {0.,0.};
  vec_t fRn = {0.,0.};
  vec_t fRp = {0.,0.};

  pressiodemoapps::weno3(fLn, fLp, fRn, fRp, q, 0, 2, 6, 8, 4, 2);

  for (int k=0; k<2; ++k)
  {
    const auto & goldF = (k==0) ? goldF0 : goldF1;

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


int main(int argc, char *argv[])
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

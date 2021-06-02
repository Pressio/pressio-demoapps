
#include "weno.hpp"
#include <array>
#include <iostream>
#include <iomanip>

int main(int argc, char *argv[])
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
    {0.498169371091894, 0.498239224509689,
     1.50249211669906, 1.53308488724067};

  for (int i=0; i<4; ++i){
    const auto err = std::abs(goldF[i] - f[i]);
    if (err > 1e-13){
      std::puts("FAILED");
      return 0;
    }
  }
  std::puts("PASS");

  return 0;
}

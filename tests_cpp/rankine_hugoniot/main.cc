
#include <array>
#include <iomanip>
#include "pressiodemoapps/euler1d.hpp"

int main()
{
  using scalar_t = double;
  using vec_t = std::array<scalar_t, 3>;
  scalar_t gamma = 1.4;
  scalar_t machShock = 10.0;
  vec_t primPreShock{gamma, 0.0, 1.0};
  vec_t primPostShock{0.0, 0.0, 0.0};

  pressiodemoapps::ee::computePostShockConditions(primPostShock,
                                                  primPreShock,
                                                  machShock,
                                                  gamma);

  for (auto & it : primPostShock){
    std::cout << std::setprecision(15) << it << " \n";
  }

  const std::array<scalar_t, 3> goldPrimPostShock{8.0, -8.25, 116.5};

  for (int i=0; i<3; ++i){
    const auto err = std::abs(goldPrimPostShock[i] - primPostShock[i]);
    if (err > 1e-13){
      std::puts("FAILED");
      return 0;
    }
  }
  std::puts("PASS");

  return 0;
}

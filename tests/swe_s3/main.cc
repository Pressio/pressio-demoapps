
#include "swe.hpp"
#include <iomanip>

int main(int argc, char *argv[])
{
  using scalar_t = double;
  using app_t       = pressiodemoapps::Swe;
  using app_state_t = typename app_t::state_type;
  using app_rhs_t   = typename app_t::velocity_type;

  std::array<scalar_t,3> params = {0.25, 0.5, 0.1};
  app_t appObj(".", params);

  auto fomVelo = appObj.createVelocity();
  app_state_t fomState(fomVelo.size());
  for (int i=0; i<fomState.size(); ++i){
    fomState(i) = (double) (i+1);
  }

  appObj.velocity(fomState, 0., fomVelo);
  //std::cout << std::setprecision(14) << fomVelo << std::endl;

  app_state_t goldVelo(fomVelo.size());
  goldVelo <<
    440.19802377938, 460.99843957979, 739.39885538021,
    356.61259942078, 346.77688513507, 696.24117084935,
    369.02281312333, 349.10531006335, 742.18780700338,
    344.18175264043, 355.57260526606, 749.6634578917,
    48.461276134967, 110.76718522588, -98.526905683215,
    -44.340216773232, -89.634490920138, -232.82876506704,
    -44.468414313835, -98.89192709138, -269.51543986893,
    -98.573684863215, -47.93991491799, -359.80614497276,
    66.592185379923, 162.47916453749, -230.43385630494,
    -46.319478775087, -127.15399490412, -379.08851103315,
    -46.880336994272, -136.45022889228, -415.42012079029,
    -114.93204009409, -27.252662304976, -519.27328451586,
    -183.7999939098, -54.746515648928, -31.693037388059,
    -310.94451610532, -427.59411653631, -116.54371696731,
    -321.93128757954, -446.86716215377, -90.403036728008,
    -412.87868107042, -288.36857643775, -144.75847180507;

  if (!goldVelo.isApprox(fomVelo)){
    std::puts("FAIL");
  }
  else{
    std::puts("PASS");
  }

  return 0;
}

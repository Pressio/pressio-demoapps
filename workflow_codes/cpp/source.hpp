
#ifndef TWODIM_REAC_Diff_SOURCE_HPP_
#define TWODIM_REAC_Diff_SOURCE_HPP_

#include <type_traits>

auto mySource = [](const double & x,
                   const double & y,
                   const double & t,
                   double & result)
{
  result = std::sin(M_PI*x*(y-0.2)) * 0.1*std::sin(4.*M_PI*y*x);

  // result = 0.25*std::sin(4.*M_PI*(x-0.2)*(y-0.5)) 
  //       + 0.001*std::exp(2.*(x-0.8)*y)
  //       + 0.001*std::exp(2.*(x-0.8)*(y-0.2));
};

#endif
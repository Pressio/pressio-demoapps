
#ifndef PRESSIODEMOAPPS_TESTS_OBSERVER_HPP_
#define PRESSIODEMOAPPS_TESTS_OBSERVER_HPP_

#include "pressio_ode_explicit.hpp"

template <typename state_t>
class FomObserver
{
public:
  FomObserver(const std::string & f0, int freq)
    : myfile0_(f0,  std::ios::out | std::ios::binary),
      sampleFreq_(freq){}

  ~FomObserver(){
    myfile0_.close();
  }

  template<typename time_t>
  void operator()(const size_t step,
  		  const time_t t,
  		  const state_t & y)
  {
    if (step % sampleFreq_ == 0){
      const auto & ydata = *y.data();
      const std::size_t ext = y.extent(0)*sizeof(double);
      myfile0_.write(reinterpret_cast<const char*>(&ydata(0)), ext);
    }
  }

private:
  std::ofstream myfile0_;
  int sampleFreq_ = {};
};

template<class T>
using RomObserver = FomObserver<T>;

#endif

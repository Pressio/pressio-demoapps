
#ifndef PRESSIODEMOAPPS_TESTS_OBSERVER_HPP_
#define PRESSIODEMOAPPS_TESTS_OBSERVER_HPP_

template <typename StateType>
class FomObserver
{
public:
  FomObserver(const std::string & f0, int freq)
    : myfile0_(f0,  std::ios::out | std::ios::binary),
      sampleFreq_(freq){}

  ~FomObserver(){
    myfile0_.close();
  }

  FomObserver() = default;
  FomObserver & operator=(FomObserver &) = default;
  FomObserver & operator=(FomObserver &&) = default;

  template<typename TimeType>
  void operator()(const pressio::ode::StepCount stepIn,
  		  const TimeType /*timein*/,
  		  const StateType & state)
  {
    const auto step = stepIn.get();
    if (step % sampleFreq_ == 0){
      const std::size_t ext = state.size()*sizeof(double);
      myfile0_.write(reinterpret_cast<const char*>(&state(0)), ext);
    }
  }

private:
  std::ofstream myfile0_;
  int sampleFreq_ = {};
};

#endif

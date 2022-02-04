
#ifndef OBSERVER_HPP_
#define OBSERVER_HPP_

template <typename ObservableType>
class StateObserver
{
public:
  StateObserver(const std::string & f0, int freq)
    : myfile_(f0,  std::ios::out | std::ios::binary),
      sampleFreq_(freq){}

  ~StateObserver(){
    myfile_.close();
  }

  template<typename time_t>
  void operator()(const size_t step,
  		  const time_t t,
  		  const ObservableType & y)
  {
    if (step % sampleFreq_ == 0){
      const std::size_t ext = y.size()*sizeof(typename ObservableType::Scalar);
      myfile_.write(reinterpret_cast<const char*>(&y(0)), ext);
    }
  }

private:
  std::ofstream myfile_;
  const int sampleFreq_ = {};
};

template <typename T = void>
class VelocityObserver
{
public:
  VelocityObserver(const std::string & f0, int freq)
    : myfile_(f0,  std::ios::out | std::ios::binary),
      sampleFreq_(freq){}

  ~VelocityObserver(){
    myfile_.close();
  }

  template<class time_t, class scalar_type>
  void operator()(const size_t step,
  		  const time_t t,
  		  const Eigen::Matrix<scalar_type, -1, 1> & f)
  {
    if (step % sampleFreq_ == 0){
      const std::size_t ext = f.size()*sizeof(scalar_type);
      myfile_.write(reinterpret_cast<const char*>(&f(0)), ext);
    }
  }

private:
  std::ofstream myfile_;
  const int sampleFreq_ = {};
};

#endif

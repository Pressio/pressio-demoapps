
#ifndef PRESSIO_DEMOAPPS_NOOP_FUNCTOR_HPP_
#define PRESSIO_DEMOAPPS_NOOP_FUNCTOR_HPP_

namespace pressiodemoapps{
namespace impl{

template <class ReturnType>
class NoOperation{
public:
  explicit NoOperation(ReturnType return_value = ReturnType())
    : value_(return_value) { }

  template<class ...Types>
  ReturnType operator () (Types && ... /*args*/) const noexcept { return value_; }

  template<class ...Types>
  ReturnType operator () (Types && ... /*args*/) noexcept { return value_; }

private:
  ReturnType value_;
};

//! Specialized noop functor which returns a void.
template <>
class NoOperation<void>{
public:
  template<class ...Types>
  void operator () (Types && ... /*args*/) const noexcept{}

  template<class ...Types>
  void operator () (Types && ... /*args*/) noexcept{}
};

}}

#endif

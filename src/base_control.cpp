#include "base_control.h"

namespace nonlinear_control {

template <typename T>
BaseController<T>::BaseController() {}

template <typename T>
BaseController<T>::~BaseController() {}

template <typename T>
void BaseController<T>::init() {}

template <typename T>
void BaseController<T>::run() {}

template class BaseController<float>;
template class BaseController<double>;

}  // namespace nonlinear_control

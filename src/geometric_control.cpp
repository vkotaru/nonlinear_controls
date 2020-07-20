#include "geometric_control.h"

namespace nonlinear_control {

template <typename T>
GeometricController<T>::GeometricController(/* args */) : BaseController<T>() {
}

template <typename T>
GeometricController<T>::~GeometricController() {
}

template <typename T>
void GeometricController<T>::init() {
}

template <typename T>
void GeometricController<T>::run() {}

/*
 * Explicitly instantiate the template, and its member definitions
 * https://stackoverflow.com/questions/8752837/undefined-reference-to-template-class-constructor
*/
template class GeometricController<float>;
template class GeometricController<double>;

}  // namespace nonlinear_control

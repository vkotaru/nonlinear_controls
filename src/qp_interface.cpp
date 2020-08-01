#include "qp_interface.h" 

namespace nonlinear_control {

template <typename T>
QPInterface<T>::QPInterface() {

}

template <typename T>
QPInterface<T>::~QPInterface() {

}

template <typename T>
void QPInterface<T>::setup() {
    
}

template <typename T>
void QPInterface<T>::solve() {
    
}

template class QPInterface<float>;
template class QPInterface<double>;

}


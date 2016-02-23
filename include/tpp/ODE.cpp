// ODE.cpp
// Implementation for abstract base class
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#include <SDE.hpp>


// Constructor
template <typename T>
ODE<T>::ODE( const int d )
  : dim( d )
{
  // empty
}

// Get the dimensions of the equation
template <typename T>
int ODE<T>::getDim() const
{
  return dim;
}


// explicit class template instatiation
template class ODE<float>;
template class ODE<double>;

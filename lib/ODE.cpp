// ODE.cpp
// Implementation for abstract base class
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#include <SDE.hpp>


// Constructor
ODE::ODE( const int d )
  : dim( d )
{
  // empty
}

// Get the dimensions of the equation
int ODE::getDim() const
{
  return dim;
}

// LangevinEquation.cpp
// Implementation for abstract base class
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#include <LangevinEquation.h>

// Constructor
LangevinEquation::LangevinEquation(int d)
  : dim(d)
{
  // empty
}

// Get the dimensions of the equation
int LangevinEquation::getDim()
{
  return dim;
}

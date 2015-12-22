// LangevinEquation.cpp
// Implementation for abstract base class
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#include <LangevinEquation.hpp>

// Langevin Equation object is required to have a fixed dimension,
// specify 'd' when constructing an instance.

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

// Langevin Equation might only have drift 
void LangevinEquation::computeDiffusion( matrix_f& out,
                                         array_f& /*no input*/,
                                         float /*no time*/ )
{
  std::fill( out.data(), out.data()+out.num_elements(), 0 );
}

// Langevin equation can be defined without specifying derivatives
void LangevinEquation::computeDiffusionDerivatives( array3_f &out, array_f&,
						    float )
{
  std::fill( out.data(), out.data()+out.num_elements(), 0 );
}

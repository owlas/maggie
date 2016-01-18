// SDE.cpp
// Implementation for abstract base class
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#include <SDE.hpp>

// Langevin Equation object is required to have a fixed dimension,
// specify 'd' when constructing an instance.

// Constructor
SDE::SDE( const int d, const int m )
    : ODE( d )
    , wDim( m )
{
    // empty
}

// Get the dimensions of the equation
int SDE::getWDim() const
{
    return wDim;
}

// Langevin equation can be defined without specifying derivatives
void SDE::computeDiffusionDerivatives( array3_f &out,
                                       const array_f&,
                                       const float ) const
{
    std::fill( out.data(), out.data()+out.num_elements(), 0 );
}

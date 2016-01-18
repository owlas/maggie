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
template <typename T>
SDE<T>::SDE( const int d, const int m )
    : ODE<T>( d )
    , wDim( m )
{
    // empty
}

// Get the dimensions of the equation
template <typename T>
int SDE<T>::getWDim() const
{
    return wDim;
}

// Langevin equation can be defined without specifying derivatives
template <typename T>
void SDE<T>::computeDiffusionDerivatives( array3<T> &out,
                                          const array<T>&,
                                          const T ) const
{
    std::fill( out.data(), out.data()+out.num_elements(), 0 );
}

// Explicit instatiation of template classes
template class SDE<float>;
template class SDE<double>;

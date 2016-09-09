// SDE.cpp
// Implementation for abstract base class
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#include <SDE.hpp>

// Langevin Equation object is required to have a fixed dimension,
// specify 'd' when constructing an instance.

// Langevin equation can be defined without specifying derivatives
template <size_t DIM, size_t WDIM>
void SDE<DIM,WDIM>::computeDiffusionDerivatives( typename SDE<DIM,WDIM>::array3 &out,
                                                 const typename SDE<DIM,WDIM>::array&,
                                                 const double ) const
{
    for( auto &dim1 : out )
        for( auto &dim2 : dim1 )
            for( auto &val : dim2 )
                val = 0;
}

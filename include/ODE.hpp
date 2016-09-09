// ODE.h
// Header for abstract base class
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef ODE_H
#define ODE_H

#include <cstddef>
#include <array>
template <size_t DIM>
using darray = std::array<double,DIM>;

template<size_t DIM>
class ODE
{
public:
    ODE() {};

    // Useful aliases
    using array = darray<DIM>;
    static const size_t dim = DIM;

    // Return derivative from a given state vector
    virtual void computeDrift( array& out, const array& in,
                               const double t ) const = 0;
};

#endif

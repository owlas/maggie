// ODE.h
// Header for abstract base class
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef ODE_H
#define ODE_H

#include <algorithm>
#include <boost/multi_array.hpp>
template <typename T> using array = boost::multi_array<T,1>;
template <typename T> using matrix = boost::multi_array<T,2>;

template <typename T>
class ODE
{
public:
    ODE( const int dim );

    int getDim() const; // get dimensions of equation

    // Return derivative from a given state vector
    virtual void computeDrift( array<T>& out, const array<T>& in,
                               const T t ) const = 0;
private:
    const int dim;
};

#include <tpp/ODE.cpp>
#endif

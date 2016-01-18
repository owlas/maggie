// SDE.h
// Header for abstract base class
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef SDE_H
#define SDE_H

#include <algorithm>
#include <ODE.hpp>
#include <boost/multi_array.hpp>
template <typename T> using array = boost::multi_array<T,1>;
template <typename T> using matrix = boost::multi_array<T,2>;
template <typename T> using array3 = boost::multi_array<T,3>;

template <typename T>
class SDE : public ODE<T>
{
public:
    SDE( const int dim, const int wDim );

    int getWDim() const; // get dimensions of equation

    // Return Langevin components from a given state vector
    virtual void computeDiffusion( matrix<T>& out, const array<T>& in,
                                   const T t ) const = 0;

    // Langevin equations can specify derivative terms for the diffusion matrix
    // BB[i][j][k] = PD of B[i][j] w.r.t. state[k]
    //
    virtual void computeDiffusionDerivatives( array3<T> &out, const array<T> &in,
                                              const T t ) const;

private:
    const int wDim;
};

#include <tpp/SDE.cpp>
#endif

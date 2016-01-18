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
using array_f = boost::multi_array<float,1>;
using matrix_f = boost::multi_array<float,2>;
using array3_f = boost::multi_array<float,3>;

class SDE : public ODE
{
public:
    SDE( const int dim, const int wDim );

    int getWDim() const; // get dimensions of equation

    // Return Langevin components from a given state vector
    virtual void computeDiffusion( matrix_f& out, const array_f& in,
                                   const float t ) const = 0;

    // Langevin equations can specify derivative terms for the diffusion matrix
    // BB[i][j][k] = PD of B[i][j] w.r.t. state[k]
    //
    virtual void computeDiffusionDerivatives( array3_f &out, const array_f &in,
                                              const float t ) const;

private:
    const int wDim;
};

#endif

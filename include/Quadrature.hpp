// Quadrature.hpp - class to integrate using quadrature rule
//
// Based on numerical recipes version
//
// Oliver W. Laslett
// O.Laslett@soton.ac.uk
//
#ifndef QUAD_H
#define QUAD_H

#include<functional>
#include<stdexcept>
#include<cmath>
#include<boost/multi_array.hpp>
using array_f=boost::multi_array<float,1>;

namespace Quad 
{

  class Quadrature
  {
  public:
    // constructor
    // specifiy a function to integrate
    Quadrature( std::function<float( float )> func, const float start,
                const float end );

    // second constructor
    // specify a vector of values to integrate
    Quadrature( array_f vec );

    // Returns the value of the integral at the nth stage of
    // refinement. 
    float next();

    // Returns n
    int getN() const;

    // Integrates with increasing refinement until a limit is reached
    float qTrap( const float eps=1e-7 );
  private:
    std::function<float( float )> f; // function to integrate
    const float a; // upper limit
    const float b; // lower limit
    float s; // current integral value
    int n; // current integral step
  };

  // A single quadrature algorithm for discrete data
  // integrates from x = start to x = end at intervals of h
  // and vec = y( x ) for start, start+h, start+2h, ...
  float trapVec( array_f vec, const float h );
}
#endif


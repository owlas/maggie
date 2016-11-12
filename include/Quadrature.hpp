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
template <typename T> using array = boost::multi_array<T,1>;

namespace Quad
{

    template <typename T>
  class Quadrature
  {
  public:
    // constructor
    // specifiy a function to integrate
    Quadrature( std::function<T( T )> func, const T start,
                const T end );

    // second constructor
    // specify a vector of values to integrate
      Quadrature( array<T> vec );

    // Returns the value of the integral at the nth stage of
    // refinement.
    T next();

    // Returns n
    int getN() const;

    // Integrates with increasing refinement until a limit is reached
    T qTrap( const T eps=1e-7 );
  private:
    std::function<T ( T )> f; // function to integrate
    const T a; // upper limit
    const T b; // lower limit
    T s; // current integral value
    int n; // current integral step
  };

  // A single quadrature algorithm for discrete data
  // integrates from x = start to x = end at intervals of h
  // and vec = y( x ) for start, start+h, start+2h, ...
    template <class Container>
    double trapVec( const Container vec, const double h );

    template <class Container>
    double trapVec( const Container yvec, const Container xvec, const int Npoints );

    template <class Container>
    double trapVec( const Container yvec, const Container xvec );
}

#include<tpp/Quadrature.cpp>
#endif

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
typedef boost::multi_array<float,1> array_f;
typedef boost::multi_array<float,2> matrix_f;

class ODE
{
 public:
  ODE( const int dim );
  
  int getDim() const; // get dimensions of equation
  
  // Return derivative from a given state vector
  virtual void computeDrift( array_f& out, const array_f& in,
			     const float t ) const = 0;
 private:
  const int dim;
};
#endif
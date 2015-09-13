// LangevinEquation.h
// Header for abstract base class
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef LANGEVIN_H
#define LANGEVIN_H

#include <boost/multi_array.hpp>
typedef boost::multi_array<float,1> array_f;
typedef boost::multi_array<float,2> matrix_f;

class LangevinEquation
{
 public:
  LangevinEquation( int dim );
  
  int getDim(); // get dimensions of equation
  
  // Return Langevin components from a given state vector
  virtual void computeDrift( array_f&, array_f& ) const = 0; // returns the drift vector
  virtual void computeDiffusion( matrix_f&, array_f& ) = 0; // returns the diffusion marix

 private:
  const int dim;
};
#endif

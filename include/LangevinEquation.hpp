// LangevinEquation.h
// Header for abstract base class
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef LANGEVIN_H
#define LANGEVIN_H

#include <algorithm>
#include <boost/multi_array.hpp>
typedef boost::multi_array<float,1> array_f;
typedef boost::multi_array<float,2> matrix_f;
typedef boost::multi_array<float,3> array3_f;

class LangevinEquation
{
 public:
  LangevinEquation( int dim );
  
  int getDim(); // get dimensions of equation
  
  // Return Langevin components from a given state vector
  virtual void computeDrift( array_f& out, array_f& in, float t ) = 0;
  virtual void computeDiffusion( matrix_f& out, array_f& in, float t ); 

  // Langevin equations can specify derivative terms for the diffusion matrix
  // BB[i][j][k] = PD of B[i][j] w.r.t. state[k]
  // 
  virtual void computeDiffusionDerivatives( array3_f &out, array_f &in,
					    float t );

 private:
  const int dim;
};
#endif

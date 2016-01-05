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
  LangevinEquation( const int dim );
  
  int getDim() const; // get dimensions of equation
  
  // Return Langevin components from a given state vector
  virtual void computeDrift( array_f& out, const array_f& in,
			     const float t ) const = 0;
  virtual void computeDiffusion( matrix_f& out, const array_f& in,
				 const float t ) const; 

  // Langevin equations can specify derivative terms for the diffusion matrix
  // BB[i][j][k] = PD of B[i][j] w.r.t. state[k]
  // 
  virtual void computeDiffusionDerivatives( array3_f &out, const array_f &in,
					    const float t ) const;

 private:
  const int dim;
};
#endif

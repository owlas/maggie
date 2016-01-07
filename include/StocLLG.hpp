// StocLLG.h
// Header for the stochastic LLG Langevin equation
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef STOCLLG_H
#define STOCLLG_H

#include <SDE.hpp>

#include <boost/multi_array.hpp>
typedef boost::multi_array<float,1> array_f;
typedef boost::multi_array<float,2> matrix_f;

class StocLLG : public SDE
{
 public:
        StocLLG( const float s, const float a,
		 const float hx, const float hy, const float hz );

  // returns the drift vector
  virtual void computeDrift( array_f&, const array_f&, const float ) const;
  // returns the diffusion marix
  virtual void computeDiffusion( matrix_f&, const array_f&,
				 const float ) const;
  // returns a vector of derivative sums for Taylor
  virtual void computeDiffusionDerivatives( array3_f &out, const array_f &in,
					    const float ) const;

  void setReducedHeff( const float, const float, const float );
  array_f getReducedHeff() const;

  void setSigma( const float );
  float getSigma() const;
  
  void setAlpha( const float );
  float getAlpha() const;

 private:
  array_f h;
  float sigma;
  float alpha;
};
#endif

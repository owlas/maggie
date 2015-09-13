// StocLLG.h
// Header for the stochastic LLG Langevin equation
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef STOCLLG_H
#define STOCLLG_H

#include <LangevinEquation.hpp>

#include <boost/multi_array.hpp>
typedef boost::multi_array<float,1> array_f;
typedef boost::multi_array<float,2> matrix_f;

class StocLLG : public LangevinEquation
{
 public:
        StocLLG( float s, float a, float M, float hx, float hy, float hz );

  virtual void computeDrift( array_f&, array_f& ) const;  // returns the drift vector
  virtual void computeDiffusion( matrix_f&, array_f& ); // returns the diffusion marix

  void setReducedHeff( float, float, float );
  array_f& getReducedHeff();

  void setSigma( float );
  float getSigma();
  
  void setAlpha( float );
  float getAlpha();

  void setMs( float );
  float getMs();

 private:
  array_f h;
  float sigma;
  float alpha;
  float Ms;
};
#endif

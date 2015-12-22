// SDE.hpp -- header file containing implementations of various stochastic
// differential equations as Langevin Equation objects.
// 
// O.Laslett@soton.ac.uk <O.Laslett@soton.ac.uk>
// 2015 
// 
#ifndef SDE_H
#define SDE_H

#include <LangevinEquation.hpp>
#include <boost/multi_array.hpp>
using array_f = boost::multi_array<float,1>;
using matrix_f = boost::multi_array<float,2>;

// Ornstein-Uhlenbeck process
class OH : public LangevinEquation
{
public:
  // constructor
  OH( float theta, float mu, float sigma );
  virtual void computeDrift( array_f& out, array_f& in, float );
  virtual void computeDiffusion( matrix_f& out, array_f& in, float );

  float getTheta() const;
  float getMu() const;
  float getSigma() const;
  
private:
  float theta, mu, sigma;
};

// Wiener process
class Wiener : public LangevinEquation
{
public:
  Wiener();
  virtual void computeDrift( array_f& out, array_f&, float );
  virtual void computeDiffusion( matrix_f& out, array_f&, float );
};

// Simple deterministic ODE with constant drift
class ODEConstantDrift : public LangevinEquation
{
public:
  ODEConstantDrift( float );
  virtual void computeDrift( array_f& out, array_f&, float );

private:
  float a;
};

// Stochastic differential equation with constant multiplicative noise
// dX = aX dt + bX dW
class MultiplicativeConstantNoise : LangevinEquation
{
public:
  MultiplicativeConstantNoise( float a, float b );
  virtual void computeDrift( array_f& out, array_f& in, float );
  virtual void computeDiffusion( matrix_f& out, array_f& in, float );
  virtual void computeDiffusionDerivatives( array3_f& out, array_f& in,
					    float );

private:
  float a, b;
};
#endif

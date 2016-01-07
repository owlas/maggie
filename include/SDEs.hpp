// SDES.hpp -- header file containing implementations of various stochastic
// differential equations as Langevin Equation objects.
// 
// O.Laslett@soton.ac.uk <O.Laslett@soton.ac.uk>
// 2015 
// 
#ifndef SDES_H
#define SDES_H

#include <SDE.hpp>
#include <boost/multi_array.hpp>
using array_f = boost::multi_array<float,1>;
using matrix_f = boost::multi_array<float,2>;

// Ornstein-Uhlenbeck process
class OH : public SDE
{
public:
  // constructor
  OH( const float theta, const float mu, const float sigma );
  virtual void computeDrift( array_f& out, const array_f& in,
			     const float ) const;
  virtual void computeDiffusion( matrix_f& out, const array_f& in,
				 const float ) const;

  float getTheta() const;
  float getMu() const;
  float getSigma() const;
  
private:
  const float theta, mu, sigma;
};

// Wiener process
class Wiener : public SDE
{
public:
  Wiener();
  virtual void computeDrift( array_f& out, const array_f&,
			     const float ) const;
  virtual void computeDiffusion( matrix_f& out, const array_f&,
				 const float ) const;
};

// Simple deterministic ODE with constant drift
class ODEConstantDrift : public ODE
{
public:
  ODEConstantDrift( const float );
  virtual void computeDrift( array_f& out, const array_f&,
			     const float ) const;

private:
  const float a;
};


// Scalar SDE with constant drift and additive noise
// dx = ax * dt + b * dW
class SDE_AXpB : public SDE
{
public:
  SDE_AXpB( const float a, const float b );
  virtual void computeDrift( array_f& out, const array_f&,
			     const float ) const;
  virtual void computeDiffusion( matrix_f& out, const array_f&,
				 const float ) const;

private:
  const float a, b;
};

// Scalar SDE with constant drift and multiplicative noise
// dx = ax * dt + bx * dW
class SDE_AXpBX : public SDE
{
public:
  SDE_AXpBX( const float a, const float b );
  virtual void computeDrift( array_f& out, const array_f&,
			     const float ) const;
  virtual void computeDiffusion( matrix_f& out, const array_f&,
				 const float ) const;

private:
  const float a, b;
};

// Stochastic differential equation with constant multiplicative noise
// dX = aX dt + bX dW
class MultiplicativeConstantNoise : SDE
{
public:
  MultiplicativeConstantNoise( const float a, const float b );
  virtual void computeDrift( array_f& out, const array_f& in,
			     const float ) const;
  virtual void computeDiffusion( matrix_f& out, const array_f& in,
				 const float ) const;
  virtual void computeDiffusionDerivatives( array3_f& out, const array_f& in,
					    const float ) const;

private:
  const float a, b;
};
#endif

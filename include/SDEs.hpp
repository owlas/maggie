// SDES.hpp -- header file containing implementations of various stochastic
// differential equations as Langevin Equation objects.
//
// O.Laslett@soton.ac.uk <O.Laslett@soton.ac.uk>
// 2015
//
#ifndef SDES_H
#define SDES_H

#include <SDE.hpp>

// Ornstein-Uhlenbeck process
class OH : public SDE<1,1>
{
public:
    // constructor
    OH( const double theta, const double mu, const double sigma );
    void computeDrift( OH::array& out, const OH::array& in,
                       const double ) const;
    void computeDiffusion( OH::matrix& out, const OH::array& in,
                           const double ) const;

    double getTheta() const;
    double getMu() const;
    double getSigma() const;

private:
    const double theta, mu, sigma;
};

// Wiener process
class Wiener : public SDE<1,1>
{
public:
    Wiener();
    void computeDrift( Wiener::array& out, const Wiener::array&,
                       const double ) const;
    void computeDiffusion( Wiener::matrix& out, const Wiener::array&,
                           const double ) const;
};

// Simple deterministic ODE with constant drift
class ODEConstantDrift : public ODE<1>
{
public:
    ODEConstantDrift( const double );
    void computeDrift( ODEConstantDrift::array& out,
                       const ODEConstantDrift::array&,
                       const double ) const;

private:
  const double a;
};


// Scalar SDE with constant drift and additive noise
// dx = ax * dt + b * dW
class SDE_AXpB : public SDE<1,1>
{
public:
    SDE_AXpB( const double a, const double b );
    void computeDrift( SDE<1,1>::array& out, const SDE<1,1>::array&,
                       const double ) const;
    void computeDiffusion( SDE<1,1>::matrix& out, const SDE<1,1>::array&,
                           const double ) const;

private:
    const double a, b;
};

// Scalar SDE with constant drift and multiplicative noise
// dx = ax * dt + bx * dW
class SDE_AXpBX : public SDE<1,1>
{
public:
  SDE_AXpBX( const double a, const double b );
  void computeDrift( SDE<1,1>::array& out, const SDE<1,1>::array&,
                     const double ) const;
  void computeDiffusion( SDE<1,1>::matrix& out, const SDE<1,1>::array&,
                         const double ) const;

private:
    const double a, b;
};

// Stochastic differential equation with constant multiplicative noise
// dX = aX dt + bX dW
class MultiplicativeConstantNoise : SDE<1,1>
{
public:
    MultiplicativeConstantNoise( const double a, const double b );
    void computeDrift( SDE<1,1>::array& out, const SDE<1,1>::array& in,
                       const double ) const;
    void computeDiffusion( SDE<1,1>::matrix& out, const SDE<1,1>::array& in,
                           const double ) const;
    void computeDiffusionDerivatives( SDE<1,1>::array3& out,
                                      const SDE<1,1>::array& in,
                                      const double ) const;

private:
    const double a, b;
};
#endif

// TwoStateMasterEquation.hpp
// Implementation of the differential equation for a two state master
// equation.
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef TWOSTATEME_H
#define TWOSTATEME_H

#include<LangevinEquation.hpp>

class TwoStateMasterEquation : public LangevinEquation
{
public:

  // constructor
  TwoStateMasterEquation( float rate1, float rate2 );

  // Compute the drift vecctor for the current state
  virtual void computeDrift( array_f& out, array_f& in );

  // set the transition rates
  void setRates( float rate1, float rate2 );
  // get rate one and two
  float getRate1();
  float getRate2();

private:
  float w1;
  float w2;
};

#endif

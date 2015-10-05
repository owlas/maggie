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
#include<functional>

class TwoStateMasterEquation : public LangevinEquation
{
public:

  // constructor for fixed rates
  TwoStateMasterEquation( float rate1, float rate2 );

  // constructor for variable rates
  TwoStateMasterEquation( std::function<float( float )> rateFunc1,
                          std::function<float( float )> rateFunc2 );

  // Compute the drift vecctor for the current state
  virtual void computeDrift( array_f& out, array_f& in, float t );

  // get rate one and two at a set time
  float getRate1( float t=0);
  float getRate2( float t=0);

private:
  std::function<float( float )> w1;
  std::function<float( float )> w2;
};

#endif

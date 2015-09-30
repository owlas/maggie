// TwoStateMasterEquation.cpp
//
// Implementation of the two dimensional ordinary differential
// equation governing the probability evolution of the two state
// master equation.
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//

#include<TwoStateMasterEquation.hpp>
#include<stdexcept>
using std::invalid_argument;

// constructor
TwoStateMasterEquation::TwoStateMasterEquation( float rate1,
						float rate2 )
  : LangevinEquation( 2 )
{
  setRates( rate1, rate2 );
}

// setter and getters for the rates
void TwoStateMasterEquation::setRates( float rate1, float rate2 )
{
  w1 = rate1;
  w2 = rate2;
}

float TwoStateMasterEquation::getRate1() { return w1; }

float TwoStateMasterEquation::getRate2() { return w2; }

// Compute the two-state master equation
void TwoStateMasterEquation::computeDrift( array_f& out, array_f& in, float )
{
  if( in[0] + in[1] == 1.0 )
    throw invalid_argument( "Error: Input vector should be a"
			    " probability vector => sum of elements"
			    " must be equal to unity" );
  out[0] = ( 1 - in[0] )*w2 - in[0]*w1;
  out[1] = 1- out[0];
}

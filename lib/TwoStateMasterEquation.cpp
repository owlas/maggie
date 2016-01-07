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
#include<iostream>
using std::cout;
using std::endl;
using std::invalid_argument;

// constructor
TwoStateMasterEquation::TwoStateMasterEquation( const float rate1,
						const float rate2 )
  : ODE( 2 )
  , w1( [ rate1 ]( float ){ return rate1; } )
  , w2( [ rate2 ]( float ){ return rate2; } )
{
  // empty
}

TwoStateMasterEquation::TwoStateMasterEquation
( const std::function<float( float )> rateFunc1,
  const std::function<float( float )> rateFunc2 )
  : ODE( 2 )
  , w1( rateFunc1 )
  , w2( rateFunc2 )
{
  // empty
}


// get the rates at a certain time
float TwoStateMasterEquation::getRate1( const float t ) const
{ return w1( t ); }

float TwoStateMasterEquation::getRate2( const float t ) const
{ return w2( t ); }

// Compute the two-state master equation
void TwoStateMasterEquation::computeDrift( array_f& out, const array_f& in,
                                           const float t ) const
{
  out[0] = -( w1( t )+w2( t ) )*in[0] + w2( t );
  out[1] = -( w2( t )+w2( t ) )*in[1] + w1( t );
}

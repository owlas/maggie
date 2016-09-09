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

// constructor
TwoStateMasterEquation::TwoStateMasterEquation( const double rate1,
						const double rate2 )
    : ODE<2>()
    , w1( [ rate1 ]( double ){ return rate1; } )
    , w2( [ rate2 ]( double ){ return rate2; } )
{
    // empty
}

TwoStateMasterEquation::TwoStateMasterEquation
( const std::function<double( double )> rateFunc1,
  const std::function<double( double )> rateFunc2 )
    : ODE<2>()
    , w1( rateFunc1 )
    , w2( rateFunc2 )
{
  // empty
}


// get the rates at a certain time
double TwoStateMasterEquation::getRate1( const double t ) const
{ return w1( t ); }

double TwoStateMasterEquation::getRate2( const double t ) const
{ return w2( t ); }

// Compute the two-state master equation
void TwoStateMasterEquation::computeDrift( ODE<2>::array& out,
                                           const ODE<2>::array& in,
                                           const double t ) const
{
  out[0] = -( w1( t )+w2( t ) )*in[0] + w2( t );
  out[1] = -( w2( t )+w2( t ) )*in[1] + w1( t );
}

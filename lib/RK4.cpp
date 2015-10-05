// RK4.cpp
// Implementation of the RK4 deterministic integrator
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#include<RK4.hpp>
#include<iostream>

// Constructor
RK4::RK4( LangevinEquation &le, array_f& init_state, float init_time,
	  float dt )
  : Integrator( le, init_state, init_time )
  , h( dt )
  , dim( le.getDim() )
  , k1( boost::extents[dim] )
  , k2( boost::extents[dim] )
  , k3( boost::extents[dim] )
  , k4( boost::extents[dim] )
  , tmp( boost::extents[dim] )
{}


void RK4::step()
{
  // Step 1
  int i=0;
  getLE().computeDrift( k1, state, getTime() );
  for( i=0; i<dim; i++ )
    tmp[i] = k1[i]*h/2.0 + state[i];

  // Step 2
  getLE().computeDrift( k2, tmp, getTime() + h/2.0);
  for( i=0; i<dim; i++ )
    tmp[i] = k2[i]*h/2.0 + state[i];

  // Step 3
  getLE().computeDrift( k3, tmp, getTime() + h/2.0 );
  for( i=0; i<dim; i++ )
    tmp[i] = k3[i]*h + state[i];

  // Step 4 and update state
  getLE().computeDrift( k4, tmp, getTime() + h);
  for( i=0; i<dim; i++ )
    tmp[i] = state[i]
      + ( h/6 )*( k1[i] + 2*k2[i] + 2*k3[i] + k4[i] );
  setState( tmp );

  // Update time
  setTime( getTime() + h );
}



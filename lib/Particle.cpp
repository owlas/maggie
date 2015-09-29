// Particle.cpp
//
// Implementation of the particle class, holds properties for a
// magnetic nano particle.
// 
// Oliver W. Laslett (2015)
// O.Laslett@soton.ac.uk
// 
#include<Particle.hpp>

// constructor
Particle::Particle( float g, float a, float m, float diam, float anis,
		    array_f anisAxis )
  : uea( boost::extents[3] )
  , gamma( g )
  , alpha( a )
  , ms( m )
  , k( anis )
{
  setUea( anisAxis );
  setSize( diam );
}

// setters for particle properties
void Particle::setGamma( float g ) { gamma=g; }
void Particle::setMs( float m ) { ms=m; }
void Particle::setAlpha( float a ) { alpha=a; }
void Particle::setK( float anis ) { k = anis; }

void Particle::setSize( float diam )
{
  d = diam;
  v = ( 4.0/3 )*M_PI*pow( d/2.0, 3 );
}

void Particle::setUea( array_f anisAxis )
{
  if( anisAxis.shape()[0] != 3 )
    throw invalid_argument( "Error: anisotropy axis must be of length"
			    " 3." );
  if( fabs( sqrt( pow( anisAxis[0], 2 ) + pow( anisAxis[1], 2 )
		   + pow( anisAxis[2], 2 ) ) - 1.0 ) > 1e-7 )
    throw invalid_argument( "Error: anisotropy axis must have a"
			    " magnitude of unity. Normalise vector" );

  uea = anisAxis;
}

// getters for properties
float Particle::getGamma() const { return gamma; }
float Particle::getAlpha() const { return alpha; }
float Particle::getMs() const { return ms; }
float Particle::getD() const { return d; }
float Particle::getK() const { return k; }
array_f Particle::getUea() const { return uea; }
float Particle::getV() const { return v; }

// Compute the energy barriers of the system in an aligned field. 
// happ is the reduced field intensity i.e. h/Hk
array_f
  Particle::alignedEnergyBarriers( float happ ) const
{
  array_f eBarriers( boost::extents[2] );
  if( happ>1 )
    throw invalid_argument( "currently aligned energy barriers are "
                             "only calculated for H < H_k" );
  eBarriers[0] = k*v*pow( 1-happ, 2 );
  eBarriers[1] = k*v*pow( 1+happ, 2 );

  return eBarriers;
}


// Compute the Neel-Arrhenius transition rates from energy barrier
// calculations in at a set temperature
array_f neelTransitionRates( float T, float b1, float b2 )
{
  // Assume tau_0
  float tau0 = 10e-10;

  array_f rates( boost::extents[2] );
  
  rates[0] = 1/( tau0*exp( b1/( Constants::KB*T ) ) );
  rates[1] = 1/( tau0*exp( b2/( Constants::KB*T ) ) );

  return rates;
}

  

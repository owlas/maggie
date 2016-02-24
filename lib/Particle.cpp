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
Particle::Particle( double g, double a, double m, double diam, double anis,
		    array_d anisAxis )
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
void Particle::setGamma( double g ) { gamma=g; }
void Particle::setMs( double m ) { ms=m; }
void Particle::setAlpha( double a ) { alpha=a; }
void Particle::setK( double anis ) { k = anis; }

void Particle::setSize( double diam )
{
  d = diam;
  v = ( 4.0/3 )*M_PI*pow( d/2.0, 3 );
}

void Particle::setUea( array_d anisAxis )
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
double Particle::getGamma() const { return gamma; }
double Particle::getAlpha() const { return alpha; }
double Particle::getMs() const { return ms; }
double Particle::getD() const { return d; }
double Particle::getK() const { return k; }
array_d Particle::getUea() const { return uea; }
double Particle::getV() const { return v; }

// Compute the energy barriers of the system in an aligned field.
// happ is the reduced field intensity i.e. h/Hk
array_d
  Particle::alignedEnergyBarriers( double happ ) const
{
  array_d eBarriers( boost::extents[2] );
  if( happ>1 )
    throw invalid_argument( "currently aligned energy barriers are "
                             "only calculated for H < H_k" );
  eBarriers[0] = k*v*pow( 1-happ, 2 );
  eBarriers[1] = k*v*pow( 1+happ, 2 );

  return eBarriers;
}


// Compute the Neel-Arrhenius transition rates from energy barrier
// calculations in at a set temperature
array_d
  Particle::neelTransitionRates( double T, double b1, double b2 ) const
{
  // Assume tau_0
  double tau0 = 1e-10;

  array_d rates( boost::extents[2] );

  rates[0] = 1/( tau0*exp( b1/( Constants::KB*T ) ) );
  rates[1] = 1/( tau0*exp( b2/( Constants::KB*T ) ) );

  return rates;
}

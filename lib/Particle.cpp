// Particle.cpp
//
// Implementation of the particle class, holds properties for a
// magnetic nano particle.
//
// Oliver W. Laslett (2015)
// O.Laslett@soton.ac.uk
//
#include<Particle.hpp>

using namespace maggie;

// constructor
Particle::Particle( magnetogyric _gamma, damping _alpha, double _ms,
                    diameter diam, anisotropy _k, axis anis )
    : gamma( _gamma )
    , alpha( _alpha )
    , ms( _ms )
    , k( _k )
{
  setUea( anis );
  setSize( diam );
}

// setters for particle properties
void Particle::setGamma( magnetogyric _gamma ) { gamma=_gamma; }
void Particle::setMs( double _ms ) { ms=_ms; }
void Particle::setAlpha( damping _alpha ) { alpha=_alpha; }
void Particle::setK( anisotropy _k ) { k = _k; }
void Particle::setUea( axis _uea ) { uea = _uea; }
void Particle::setSize( diameter diam )
{
  d = diam;
  v = ( 4.0/3 )*M_PI*pow( d/2.0, 3 );
}



// getters for properties
magnetogyric Particle::getGamma() const { return gamma; }
damping Particle::getAlpha() const { return alpha; }
double Particle::getMs() const { return ms; }
diameter Particle::getD() const { return d; }
anisotropy Particle::getK() const { return k; }
axis Particle::getUea() const { return uea; }
volume Particle::getV() const { return v; }

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
  Particle::neelTransitionRates( temperature T, double b1, double b2 ) const
{
  // Assume tau_0
  double tau0 = 1e-10;

  array_d rates( boost::extents[2] );

  rates[0] = 1/( tau0*exp( b1/( Constants::KB*T ) ) );
  rates[1] = 1/( tau0*exp( b2/( Constants::KB*T ) ) );

  return rates;
}


// Computes the reduced effective field of the uniaxial anisotropy
void Particle::computeAnisotropyField( field& h, anisotropy k, const moment& state ) const
{
    double mdote =
        uea[0] * state[0] +
        uea[1] * state[1] +
        uea[2] * state[2];
    h[0] = mdote * uea[0] * k;
    h[1] = mdote * uea[1] * k;
    h[2] = mdote * uea[2] * k;
}

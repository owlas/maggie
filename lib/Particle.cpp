// Particle.cpp
//
// Implementation of the particle class, holds properties for a
// magnetic nano particle.
// 
// Oliver W. Laslett (2015)
// O.Laslett@soton.ac.uk
// 
#define _USE_MATH_DEFINES
#include<cmath>
using std::sqrt;
using std::pow;
using std::fabs;
#include<Particle.hpp>
#include<stdexcept>
using std::invalid_argument;

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
float Particle::getGamma() { return gamma; }
float Particle::getAlpha() { return alpha; }
float Particle::getMs() { return ms; }
float Particle::getD() { return d; }
float Particle::getK() { return k; }
array_f Particle::getUea() { return uea; }
float Particle::getV() { return v; }

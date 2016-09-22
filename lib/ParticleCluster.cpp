// ParticleCluster.cpp
// Implementation of the particle cluster constant class for holding
// particles in groups in 3d space.
//
// Oliver W. Laslett (2015)
// O.Laslett@soton.ac.uk
//
#include "../include/ParticleCluster.hpp"
using namespace maggie;

// constructor
ParticleCluster::ParticleCluster( const Particle::vector list,
                                  const std::vector<position> locations )
  : N( int( list.size() ) )
  , dist( boost::extents[N][N][3] )
  , reduced_dist( boost::extents[N][N][3] )
{
  setParticles( list );
  setLocs( locations );

  // Compute the average K value
  k_av = 0;
  for( unsigned int i=0; i!=N; ++i )
      k_av += particles[i].getK();
  k_av /= N;

  // Compute the average volume
  v_av = 0;
  for( unsigned int i=0; i!=N; ++i )
      v_av += particles[i].getV();
  v_av /= N;

  // Compute the reduced K and V values for each particle
  for( unsigned int i=0; i!=N; ++i)
  {
      reduced_k.push_back( particles[i].getK() / k_av );
      reduced_v.push_back( particles[i].getV() / v_av );
  }

  // Normalise the distances by the average volume
  for( unsigned int i=0; i!=N; ++i )
  {
      for( unsigned int j=0; j!=N; ++j )
      {
          for( unsigned int k=0; k!=3; ++k )
          {
              reduced_dist[i][j][k] = dist[i][j][k] / std::cbrt(v_av);
          }
      }
  }
}

// Set the particle list
void ParticleCluster::setParticles( const Particle::vector p )
{
  particles = p;
}

// set the locations and check dimensions
// then compute the relative distances
void ParticleCluster::setLocs( const std::vector<position> l )
{
    if( ( l.size() != N ) )
        throw std::invalid_argument( "location vector must be size N"
                                     " where particle list is length N." );
    locs = l;

    for( array3_d::index i=0; i<N; i++ )
        for( array3_d::index j=0; j<N; j++ )
            for( array3_d::index k=0; k<3; k++ )
                dist[i][j][k] = std::abs( locs[i][k] - locs[j][k] );
}

// getters
unsigned int ParticleCluster::getNParticles() const { return N; }
std::vector<position> ParticleCluster::getLocations()  const { return locs; }
array3_d ParticleCluster::getDistances() const { return dist; }
const array3_d& ParticleCluster::getDistancesRef() const {
    return dist;
}
const array3_d& ParticleCluster::getReducedDistancesRef() const {
    return reduced_dist;
}
Particle ParticleCluster::getParticle( const int n )
    const { return particles[n]; }
std::vector<anisotropy> ParticleCluster::getReducedAnisConstants() const {
    return reduced_k;
}
std::vector<volume> ParticleCluster::getReducedVolumes() const {
    return reduced_v;
}

anisotropy ParticleCluster::getAverageAnisConstant() const { return k_av; }
volume ParticleCluster::getAverageVolume() const { return v_av; }

// Compute the stability ratio for each of the particles
std::vector<double> ParticleCluster::computeStability( temperature T ) const
{
  std::vector<double> srs;
  srs.assign( N,0 );

  for( unsigned int i=0; i<N; i++ )
    srs[i] = ( particles[i].getK() * particles[i].getV()
               / ( Constants::KB*T ) );

  return srs;
}


// Compute the energy barrier for each particle
// Very rough implementation of a continuation type method
// Requires fairly rigorous testing
// Only works with weak fields
array_d ParticleCluster::computeBarriers( double happ ) const
{
  // Maybe do a check for weak fields, ensure high energy barriers

  // Holds the state of the particle cluster
  matrix_d state( boost::extents[N][3] );

  // Possible system states
  int nStates = std::pow( N,2 );
  array3_d possibleStates( boost::extents[nStates][N][3] );

  // Find the metastable states
  for( int n=0; n<nStates; n++ )
    {

      // Initialise the state of each particle to lie on its anisotropy
      // axis, positive or negative depending on iteration
      axis uea;
      for( unsigned int i=0; i<N; i++ )
	{
	  // Each particle can be up or down, binary representation of
	  // the state number is used to determine this. 0 = down
	  // (along the negative anisotropy axis, and visa-versa
	  //
	  // Therefore shift the state number by particle number to
	  // obtain the particle parity.
	  uea = particles[i].getUea();
	  if( !( n >> i ) & 1 )
	    for( int j=0; j<3; j++ )
	      uea[j] = ( -1 )*uea[j];
	  for( int j=0; j<3; j++ )
	    state[i][j] = uea[j];
	}
      // Here do minimisation
    }
  // Just some return
  array_d somereturn( boost::extents[1] );
  somereturn[0] = happ;
  return somereturn;
}

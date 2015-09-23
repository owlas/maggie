// ParticleCluster.cpp
// Implementation of the particle cluster constant class for holding
// particles in groups in 3d space.
//
// Oliver W. Laslett (2015)
// O.Laslett@soton.ac.uk
//
#include<ParticleCluster.hpp>

// constructor
ParticleCluster::ParticleCluster( std::vector<Particle> list,
                                  matrix_f locations )
  : N( int( list.size() ) )
  , locs( boost::extents[N][3] )
  , dist( boost::extents[N][N][3] )
{
  setParticles( list );
  setLocs( locations );
}

// Set the particle list
void ParticleCluster::setParticles( std::vector<Particle> p )
{
  particles = p;
}

// set the location matrix and check dimensions
// then compute the relative distances
void ParticleCluster::setLocs( matrix_f l )
{
  if( ( l.shape()[0] != N ) or ( l.shape()[1] != 3 ) )
    throw std::invalid_argument( "location matrix must be size (Nx3)"
                                 " where particle list is length N." );
  locs = l;

  for( array3_f::index i=0; i<N; i++ )
    for( array3_f::index j=0; j<N; j++ )
      for( array3_f::index k=0; k<3; k++ )
        dist[i][j][k] = std::abs( locs[i][k] - locs[j][k] );
}

// getters
unsigned int ParticleCluster::getNParticles() const { return N; } 
matrix_f ParticleCluster::getLocations()  const { return locs; }
array3_f ParticleCluster::getDistances() const { return dist; }

// Compute the stability ratio for each of the particles
array_f ParticleCluster::computeStability( float T ) const
{
  array_f srs( boost::extents[N] );

  for( unsigned int i=0; i<N; i++ )
    srs[i] = ( particles[i].getK() * particles[i].getV() / ( KB*T ) );
  
  return srs;
}

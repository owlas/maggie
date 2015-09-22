// ParticleCluster.hpp
// Header file for a cluster of magnetic nanoparticles. Particle
// cluster is always a constant. Use this to store multiple particles
// in a group and arranged in three dimensional space.
//
// Oliver W. Laslett (2015)
// O.Laslett@soton.ac.uk
//
#ifndef PCLUSTER_H
#define PCLUSTER_H

#include<Particle.hpp>
#include<stdexcept>
#include<cmath>
#include<vector>
#include<boost/multi_array.hpp>
typedef boost::multi_array<float,2> matrix_f;
typedef boost::multi_array<float,3> array3_f;


class ParticleCluster
{
public:
  // constructor
  // Particle list could be a boost pointer array?
  // List of N particles and matrix of locations (N x 3)
  ParticleCluster( std::vector<Particle> list, matrix_f locations );

  // Set the particle list
  void setParticles( std::vector<Particle> );

  // Set the locations and compute distances
  void setLocs( matrix_f );

  // Compute the stability ratio of each particle
  std::vector<float> computeStability( float temperature ) const;

  // getters
  unsigned int getNParticles() const;
  matrix_f getLocations() const;
  array3_f getDistances() const;
  

private:
  unsigned int const N;
  const float KB = 1.3806485e-23;
  std::vector<Particle> particles;
  matrix_f locs;
  array3_f dist;
};
#endif
  

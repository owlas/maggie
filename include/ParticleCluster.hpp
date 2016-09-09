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

#include "Particle.hpp"

#include<stdexcept>
#include<cmath>
#include<vector>
#include<boost/multi_array.hpp>
using array_d = boost::multi_array<double,1>;
using matrix_d = boost::multi_array<double,2>;
using array3_d = boost::multi_array<double,3>;

#include <array>
using d3 = std::array<double, 3>;

#include "types.hpp"

class ParticleCluster
{
public:
    // constructor
    // Particle list could be a boost pointer array_d?
    // List of N particles and matrix_d of locations (N x 3)
    ParticleCluster( const Particle::vector list,
                     const std::vector<maggie::position> locations );

    // Set the particle list
    void setParticles( const Particle::vector );

    // Set the locations and compute distances
    void setLocs( const std::vector<maggie::position> );

    // Compute the stability ratio of each particle
    std::vector<double> computeStability( maggie::temperature temperature ) const;

    // Compute the energy barrier for each particle
    array_d computeBarriers( double happ ) const;

    // getters
    unsigned int getNParticles() const;
    std::vector<maggie::position> getLocations() const;

    array3_d getDistances() const;
    const array3_d& getDistancesRef() const;
    const array3_d& getReducedDistancesRef() const;

    Particle getParticle( const int n ) const;

    std::vector<maggie::anisotropy> getReducedAnisConstants() const;
    std::vector<maggie::volume> getReducedVolumes() const;

    maggie::anisotropy getAverageAnisConstant() const;
    maggie::volume getAverageVolume() const;


private:
    unsigned int const N;
    Particle::vector particles;
    std::vector<maggie::volume> reduced_v;
    std::vector<maggie::anisotropy> reduced_k;
    std::vector<maggie::position> locs;
    array3_d dist, reduced_dist;
    maggie::anisotropy k_av;
    maggie::volume v_av;
};

#endif

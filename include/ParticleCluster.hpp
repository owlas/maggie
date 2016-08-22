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


class ParticleCluster
{
public:
    // constructor
    // Particle list could be a boost pointer array_d?
    // List of N particles and matrix_d of locations (N x 3)
    ParticleCluster( const std::vector<Particle> list, const matrix_d locations );

    // Set the particle list
    void setParticles( const std::vector<Particle> );

    // Set the locations and compute distances
    void setLocs( const matrix_d );

    // Compute the stability ratio of each particle
    std::vector<double> computeStability( double temperature ) const;

    // Compute the energy barrier for each particle
    array_d computeBarriers( double happ ) const;

    // getters
    unsigned int getNParticles() const;
    matrix_d getLocations() const;

    array3_d getDistances() const;
    const array3_d& getDistancesRef() const;
    const array3_d& getReducedDistancesRef() const;

    Particle getParticle( const int n ) const;

    std::vector<double> getReducedAnisConstants() const;
    std::vector<double> getReducedVolumes() const;

    double getAverageAnisConstant() const;
    double getAverageVolume() const;


private:
    unsigned int const N;
    std::vector<Particle> particles;
    std::vector<double> reduced_v, reduced_k;
    matrix_d locs;
    array3_d dist, reduced_dist;
    double k_av, v_av;
};

#endif

// Particle.hpp
// Head file for the particle class. Represents a single magnetic
// nanoparticle with its parameters.
//
// Oliver W. Laslett (2015)
// O.Laslett@soton.ac.uk
//
#ifndef PARTICLE_H
#define PARTICLE_H

#define _USE_MATH_DEFINES
#include<cmath>
using std::sqrt;
using std::pow;
using std::fabs;
using std::exp;

#include "types.hpp"

#include<boost/multi_array.hpp>
using array_d = boost::multi_array<double,1>;

#include "Constants.hpp"

#include<stdexcept>
using std::invalid_argument;

class Particle
{
public:
    // constructor
    Particle( maggie::magnetogyric gamma, maggie::damping alpha, double Ms,
              maggie::diameter D, maggie::anisotropy K, maggie::axis uea );

    // aliases
    using vector = std::vector<Particle>;

    // getters for particle properties
    maggie::magnetogyric getGamma() const;
    maggie::damping getAlpha() const;
    double getMs() const;
    maggie::diameter getD() const;
    maggie::anisotropy getK() const;
    maggie::axis getUea() const;
    maggie::volume getV() const;

    // setters for particle properties
    void setGamma( maggie::magnetogyric );
    void setAlpha( maggie::damping );
    void setMs( double );
    void setSize( maggie::diameter );
    void setK( maggie::anisotropy );
    void setUea( maggie::axis );

    // Compute the energy barriers for the case of an applied field
    // parallel to the uniaxial anisotropy axis and with an intensity of
    // less that H_k. i.e. h = H/H_k < 1
    array_d alignedEnergyBarriers( double happ ) const;

    // Compute the Neel-Arrhenius transition rates at a given temperature
    array_d neelTransitionRates( maggie::temperature temp, double barrier1o,
                                 double barrier2 )     const;

    // Compute the effective field from the anisotropy
    void computeAnisotropyField( maggie::field&, const maggie::moment& ) const;



private:
    maggie::axis uea;
    maggie::magnetogyric gamma;
    maggie::damping alpha;
    double ms;
    maggie::diameter d;
    maggie::anisotropy k;
    maggie::volume v;
};
#endif

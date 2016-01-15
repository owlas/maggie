// Euler.hpp
// Header for the Euler numerical solver for SDEs
//
// Oliver W. Laslett <O.Laslett@soton.ac.uk>
// 2015
#ifndef EULER_H
#define EULER_H

#include<Integrator.hpp>
#include<SDE.hpp>
#include<cmath>
using std::sqrt;
#include<boost/random.hpp>
using boost::variate_generator;
using boost::mt19937;
using boost::normal_distribution;
using boost::extents;
#include<boost/multi_array.hpp>
using array_f = boost::multi_array<float,1>;
using matrix_f = boost::multi_array<float,2>;

class Euler : public Integrator<SDE>
{
public:
    Euler( const SDE& sde, const array_f& init_state, const float time,
           const float dt, mt19937& rng );

    virtual void step();

    // Manual wiener process mode
    // allow user to control internal wiener increments themselves
    void setManualWienerMode( const bool );
    void setWienerIncrements( const array_f );

private:
    const float h;
    const float dim, wDim;
    array_f dw;
    array_f tmp1;
    matrix_f tmp2;
    array_f xpred;
    normal_distribution<float> dist;
    variate_generator<mt19937 &, normal_distribution<float> > gen;
    bool manualWiener{ false };
};
#endif

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
template <typename T> using array = boost::multi_array<T,1>;
template <typename T> using matrix = boost::multi_array<T,2>;

template <typename T>
class Euler : public Integrator<SDE<T>, T>
{
public:
    Euler( const SDE<T>& sde, const array<T>& init_state, const T time,
           const T dt, mt19937& rng );

    virtual void step();

    // Manual wiener process mode
    // allow user to control internal wiener increments themselves
    void setManualWienerMode( const bool );
    void setWienerIncrements( const array<T> );

private:
    const T h;
    const T dim, wDim;
    array<T> dw;
    array<T> tmp1;
    matrix<T> tmp2;
    array<T> xpred;
    normal_distribution<T> dist;
    variate_generator<mt19937 &, normal_distribution<T> > gen;
    bool manualWiener{ false };
};

#include<tpp/Euler.cpp>
#endif

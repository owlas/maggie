// Heun.hpp
// Header for the Heun integration scheme
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef HEUN_H
#define HEUN_H

#include<Integrator.hpp>
#include<cmath>
using std::sqrt;
#include<boost/random.hpp>
using boost::variate_generator;
using boost::mt19937;
using boost::normal_distribution;
#include<boost/multi_array.hpp>
template <typename T> using array = boost::multi_array<T,1>;
template <typename T> using matrix = boost::multi_array<T,2>;

template <typename T>
class Heun : public Integrator<SDE<T>, T>
{

public:
    // Constructor
    Heun( const SDE<T> &sde, const array<T>& init_state,
          const T time, const T dt, mt19937 &rng );

    // Step the integrator once
    virtual void step();

    // Manual wiener process mode
    // allow user to control internal wiener increments themselves
    void setManualWienerMode( const bool );
    void setWienerIncrements( const array<T> );

private:
    const T h;
    const int dim, wDim;
    array<T> dw;
    array<T> xPred;
    array<T> tmp1;
    array<T> tmp1Up;
    matrix<T> tmp2;
    matrix<T> tmp2Up;
    normal_distribution<T> dist;
    variate_generator<mt19937 &, normal_distribution<T> > gen;
    bool manualWiener{ false };
};

#include <tpp/Heun.cpp>
#endif

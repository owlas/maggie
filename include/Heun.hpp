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
#include <vector>

#include "types.hpp"
using namespace maggie;

template <class EQ>
class Heun : public Integrator<EQ>
{

public:
    // Constructor
    Heun( const EQ &sde, const typename EQ::array& init_state,
          const double time, const double dt, mt19937 &rng );

    // Useful aliases
    using vector = std::vector<Heun<EQ>>;

    // Step the integrator once
    virtual void step();

    // Manual wiener process mode
    // allow user to control internal wiener increments themselves
    void setManualWienerMode( const bool );
    void setWienerIncrements( const typename EQ::array& );

private:
    const double h;
    const size_t dim, wDim;
    typename EQ::array dw;
    typename EQ::array xPred;
    typename EQ::array tmp1;
    typename EQ::array tmp1Up;
    typename EQ::matrix tmp2;
    typename EQ::matrix tmp2Up;
    normal_distribution<double> dist;
    variate_generator<mt19937 &, normal_distribution<double> > gen;
    bool manualWiener{ false };
};

#include <tpp/Heun.cpp>
#endif

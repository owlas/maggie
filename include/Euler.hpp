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

#include <vector>

template <class C>
class Euler : public Integrator<C>
{
public:
    Euler( const C& sde, const typename C::array& init_state, const double time,
           const double dt, mt19937& rng );

    virtual void step();

    // useful aliases
    using vector = std::vector<Euler<C>>;

    // Manual wiener process mode
    // allow user to control internal wiener increments themselves
    void setManualWienerMode( const bool );
    void setWienerIncrements( const typename C::array );

private:
    const double h;
    const size_t dim, wDim;
    typename C::array dw;
    typename C::array tmp1;
    typename C::matrix tmp2;
    typename C::array xpred;
    normal_distribution<double> dist;
    variate_generator<mt19937 &, normal_distribution<double> > gen;
    bool manualWiener{ false };
};

#include<tpp/Euler.cpp>
#endif

// Milstein.hpp
// Header for the Milstein numerical solver for SDEs
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef MILST_H
#define MILST_H

#include <Integrator.hpp>

#include <boost/multi_array.hpp>
using boostdmatrix = boost::multi_array<double,2>;

#include <boost/random.hpp>
using boost::variate_generator;
using boost::mt19937;
using boost::normal_distribution;
#include <StocLLG.hpp>
#include <cmath>
using std::sqrt; using std::pow;
#include <algorithm>

template <class C>
class Milstein : public Integrator<C>
{
public:
    // Constructor
    Milstein( const C &le, const typename C::array& init_state,
              const double time, const double dt,
              mt19937 &rng_1, mt19937 &rng_2 );

    // One step of the integrator
    virtual void step();

    // Fourier series approximation of double integral
    // may only be applicable for stratonovich
    void setP( const int pset );
    int getP() const;

    // Manual wiener process mode
    // allow user to control internal wiener increments themselves
    void setManualWienerMode( const bool );
    void setWienerIncrements( const typename C::array );

private:
    const double h;
    const size_t dim, wDim;
    normal_distribution<double> dist_1;
    normal_distribution<double> dist_2;
    variate_generator<mt19937&, normal_distribution<double> > gen_1;
    variate_generator<mt19937&, normal_distribution<double> > gen_2;
    typename C::array next_state;
    typename C::array dw;
    typename C::array dw2;
    typename C::array tmp1;
    typename C::matrix tmp2;
    typename C::matrix tmp22;
    typename C::array3 tmp3;
    unsigned int p; // Fourier series truncation
    typename C::array mu;
    boostdmatrix eta;
    boostdmatrix zeta;
    bool manualWiener{ false };
};

#include <tpp/Milstein.cpp>
#endif

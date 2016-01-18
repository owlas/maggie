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
template <typename T> using array = boost::multi_array<T,1>;
template <typename T> using matrix = boost::multi_array<T,2>;
template <typename T> using array3 = boost::multi_array<T,3>;

#include <boost/random.hpp>
using boost::variate_generator;
using boost::mt19937;
using boost::normal_distribution;
#include <StocLLG.hpp>
#include <cmath>
using std::sqrt; using std::pow;
#include <algorithm>

template <typename T>
class Milstein : public Integrator<SDE, T>
{
public:
    // Constructor
    Milstein( const SDE &le, const array<T>& init_state,
              const T time, const T dt,
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
    void setWienerIncrements( const array<T> );

private:
    const T h;
    const int dim, wDim;
    normal_distribution<T> dist_1;
    normal_distribution<T> dist_2;
    variate_generator<mt19937&, normal_distribution<T> > gen_1;
    variate_generator<mt19937&, normal_distribution<T> > gen_2;
    array<T> next_state;
    array<T> dw;
    array<T> dw2;
    array<T> tmp1;
    matrix<T> tmp2;
    matrix<T> tmp22;
    array3<T> tmp3;
    unsigned int p; // Fourier series truncation
    array<T> mu;
    matrix<T> eta;
    matrix<T> zeta;
    bool manualWiener{ false };
};

#include <tpp/Milstein.cpp>
#endif

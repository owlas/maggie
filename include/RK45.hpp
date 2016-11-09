// RK4.hpp
// Header for the RK45 adaptive integrator
//
// Oliver W. Laslett 2016
// O.Laslett@soton.ac.uk
//
#ifndef RK45_H
#define RK45_H
#include "Integrator.hpp"
#include <cmath>
#include <fenv.h>

template <class C>
class RK45 : public Integrator<C>
{

public:
    //Constructor
    RK45( const C &le, const typename C::array& init_state,
          const double init_time, const double eps );

    // Step the integrator once
    virtual void step();

    // update the step size
    void setStepSize( const double h );

private:
    double h;
    const unsigned dim;
    const double eps;
    typename C::array k1;
    typename C::array k2;
    typename C::array k3;
    typename C::array k4;
    typename C::array k5;
    typename C::array k6;
    typename C::array tmp;
    typename C::array tmp2;
    typename C::array relative_errs;

    // Cash-Karp parameter table
    static const double c11 = 0.2;
    static const double c21 = 3.0/40.0;
    static const double c22 = 9.0/40.0;
    static const double c31 = 3.0/10.0;
    static const double c32 = -9.0/10.0;
    static const double c33 = 6.0/5.0;
    static const double c41 = -11.0/54.0;
    static const double c42 = 2.5;
    static const double c43 = -70.0/27.0;
    static const double c44 = 35.0/27.0;
    static const double c51 = 1631.0/55296.0;
    static const double c52 = 175.0/512.0;
    static const double c53 = 575.0/13824.0;
    static const double c54 = 44275.0/110592.0;
    static const double c55 = 253.0/4096.0;

    static const double hc1 = 0.2;
    static const double hc2 = 0.3;
    static const double hc3 = 0.6;
    static const double hc4 = 1.0;
    static const double hc5 = 7.0/8.0;

    static const double x11 = 37.0/378.0;     // 5th order params
    static const double x13 = 250.0/621.0;
    static const double x14 = 125.0/594.0;
    static const double x16 = 512.0/1771.0;
    static const double x21 = 2825.0/27648.0; // 4th order params
    static const double x23 = 18575.0/48384.0;
    static const double x24 = 13525.0/55296.0;
    static const double x25 = 277.0/14336.0;
    static const double x26 = 0.25;
};

#include <tpp/RK45.cpp>
#endif

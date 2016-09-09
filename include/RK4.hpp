// RK4.hpp
// Header for the RK4 Integrator
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef RK4_H
#define RK4_H
#include<Integrator.hpp>

template <class C>
class RK4 : public Integrator<C>
{

public:
    //Constructor
    RK4( const C &le, const typename C::array& init_state,
         const double init_time, const double dt );

    // Step the integrator once
    virtual void step();

private:
    const double h;
    const unsigned dim;
    typename C::array k1;
    typename C::array k2;
    typename C::array k3;
    typename C::array k4;
    typename C::array tmp;
};

#include <tpp/RK4.cpp>
#endif

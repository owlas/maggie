// RK4.hpp
// Header for the RK4 Integrator
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef RK4_H
#define RK4_H
#include<Integrator.hpp>

template <typename T>
class RK4 : public Integrator<ODE, T>
{

public:
    //Constructor
    RK4( const ODE &le, const array<T>& init_state,
         const T init_time, const T dt );

    // Step the integrator once
    virtual void step();

private:
    const T h;
    const int dim;
    array<T> k1;
    array<T> k2;
    array<T> k3;
    array<T> k4;
    array<T> tmp;
};

#include <tpp/RK4.cpp>
#endif

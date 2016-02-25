// Integrator.hpp
// Header for the absract base class for integrators
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include<SDE.hpp>
#include<boost/multi_array.hpp>
template<typename T>
using array = boost::multi_array<T,1>;
using bidx = boost::multi_array_types::index;
#include<iostream>
using std::cout;
using std::endl;
#include<Integrator.hpp>
#include<stdexcept>
using std::invalid_argument;

template<class C, typename T>
class Integrator
{
 public:

    Integrator( const C& equation, const array<T>&,
                const T time=0 );

    // Get the current state
    const array<T>& getState() const;
    // Set the integrator state
    void setState( const array<T>& );

    // Get the current time
    T getTime() const;
    // Set the integrator time
    void setTime( const T );

    // Get the pointer to the associated Langevin Equation
    const C& getLE() const;

    // Compute an integration step
    virtual void step() = 0;

    // Reset the integrator to initial conditions
    // or a new initial condition
    void reset();
    void reset( const array<T> );

protected:
    array<T> state;

private:
    const C* lEq;
    array<T> initial_state;
    T t;
    T initial_t;
};

#include "tpp/Integrator.cpp"
#endif

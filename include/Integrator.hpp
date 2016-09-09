// Integrator.hpp
// Header for the absract base class for integrators
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include<SDE.hpp>
#include<Integrator.hpp>
#include<stdexcept>
using std::invalid_argument;

template<class EQ_C>
class Integrator
{
 public:
    Integrator( const EQ_C& equation, const typename EQ_C::array&,
                const double time=0 );

    // Get the current state
    const typename EQ_C::array& getState() const;
    // Set the integrator state
    void setState( const typename EQ_C::array& );

    // Get the current time
    double getTime() const;
    // Set the integrator time
    void setTime( const double );

    // Get the reference to the associated Langevin Equation
    const EQ_C& getLE() const;

    // Compute an integration step
    virtual void step() = 0;

    // Reset the integrator to initial conditions
    // or a new initial condition
    void reset();
    void reset( const typename EQ_C::array );

protected:
    typename EQ_C::array state;

private:
    const EQ_C* lEq;
    typename EQ_C::array initial_state;
    double t;
    double initial_t;
};

#include "tpp/Integrator.cpp"
#endif

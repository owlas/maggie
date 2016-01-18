// RK4.cpp
// Implementation of the RK4 deterministic integrator
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#include<RK4.hpp>

// Constructor
template <typename T>
RK4<T>::RK4( const ODE &le, const array<T>& init_state,
             const T init_time, const T dt )
    : Integrator<ODE, T>( le, init_state, init_time )
    , h( dt )
    , dim( le.getDim() )
    , k1( boost::extents[dim] )
    , k2( boost::extents[dim] )
    , k3( boost::extents[dim] )
    , k4( boost::extents[dim] )
    , tmp( boost::extents[dim] )
{}


template <typename T>
void RK4<T>::step()
{
    // Step 1
    int i=0;
    this->getLE().computeDrift( k1, this->state, this->getTime() );
    for( i=0; i<dim; i++ )
        tmp[i] = k1[i]*h/2.0 + this->state[i];

    // Step 2
    this->getLE().computeDrift( k2, tmp, this->getTime() + h/2.0);
    for( i=0; i<dim; i++ )
        tmp[i] = k2[i]*h/2.0 + this->state[i];

    // Step 3
    this->getLE().computeDrift( k3, tmp, this->getTime() + h/2.0 );
    for( i=0; i<dim; i++ )
        tmp[i] = k3[i]*h + this->state[i];

    // Step 4 and update state
    this->getLE().computeDrift( k4, tmp, this->getTime() + h);
    for( i=0; i<dim; i++ )
        tmp[i] = this->state[i]
            + ( h/6 )*( k1[i] + 2*k2[i] + 2*k3[i] + k4[i] );
    this->setState( tmp );

    // Update time
    this->setTime( this->getTime() + h );
}

// explicit template instantiation
template class RK4<float>;

// Euler.cpp
// Implementation for the Euler numerical solver for SDEs
//
// Oliver Laslett <O.Laslett@soton.ac.uk>
// 2015
//
#include <Euler.hpp>

// Constructor
template <typename T>
Euler<T>::Euler( const SDE& sde, const array<T>& init_state, const T time,
                 const T dt, mt19937& rng )
    : Integrator<SDE, T>( sde, init_state, time )
    , h( dt )
    , dim( sde.getDim() )
    , wDim( sde.getWDim() )
    , dw( extents[wDim] )
    , tmp1( extents[dim] )
    , tmp2( extents[dim][wDim])
    , xpred( extents[dim] )
    , dist( 0,1 )
    , gen( rng, dist )
{}

template <typename T>
void Euler<T>::step()
{
    // generate the brownian path
    if( !manualWiener )
        for( auto &i : dw )
            i = sqrt(h) * gen();

    this->getLE().computeDrift( tmp1, this->state, this->getTime() );
    this->getLE().computeDiffusion( tmp2, this->state, this->getTime() );

    for( bidx i=0; i!=dim; ++i )
    {
        xpred[i] = this->state[i] + tmp1[i]*h;
        for( bidx j=0; j!=wDim; ++j )
            xpred[i] += tmp2[i][j]*dw[j];
    }

    this->setState( xpred );
    this->setTime( this->getTime() + h );
}

// Set the manual wiener process mode
template <typename T>
void Euler<T>::setManualWienerMode( const bool s ) { manualWiener=s; }
template <typename T>
void Euler<T>::setWienerIncrements( const array_f a ) { dw=a; }

// The following types are visible
template class Euler<float>;

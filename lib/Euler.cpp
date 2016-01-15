// Euler.cpp
// Implementation for the Euler numerical solver for SDEs
//
// Oliver Laslett <O.Laslett@soton.ac.uk>
// 2015
//
#include <Euler.hpp>

// Constructor
Euler::Euler( const SDE& sde, const array_f& init_state, const float time,
              const float dt, mt19937& rng )
    : Integrator<SDE>( sde, init_state, time )
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

void Euler::step()
{
    // generate the brownian path
    if( !manualWiener )
        for( auto &i : dw )
            i = sqrt(h) * gen();

    getLE().computeDrift( tmp1, state, getTime() );
    getLE().computeDiffusion( tmp2, state, getTime() );

    for( array_f::index i=0; i!=dim; ++i )
    {
        xpred[i] = state[i] + tmp1[i]*h;
        for( array_f::index j=0; j!=wDim; ++j )
            xpred[i] += tmp2[i][j]*dw[j];
    }

    setState( xpred );
    setTime( getTime() + h );
}

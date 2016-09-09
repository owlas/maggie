// Euler.cpp
// Implementation for the Euler numerical solver for SDEs
//
// Oliver Laslett <O.Laslett@soton.ac.uk>
// 2015
//

// Constructor
template <class C>
Euler<C>::Euler( const C& sde, const typename C::array& init_state, const double time,
                 const double dt, mt19937& rng )
    : Integrator<C>( sde, init_state, time )
    , h( dt )
    , dim( C::dim )
    , wDim( C::wdim )
    , dist( 0,1 )
    , gen( rng, dist )
{}

template <class C>
void Euler<C>::step()
{
    // generate the brownian path
    if( !manualWiener )
        for( auto &i : dw )
            i = sqrt(h) * gen();

    this->getLE().computeDrift( tmp1, this->state, this->getTime() );
    this->getLE().computeDiffusion( tmp2, this->state, this->getTime() );

    for( unsigned int i=0; i!=dim; ++i )
    {
        xpred[i] = this->state[i] + tmp1[i]*h;
        for( unsigned int j=0; j!=wDim; ++j )
            xpred[i] += tmp2[i][j]*dw[j];
    }

    this->setState( xpred );
    this->setTime( this->getTime() + h );
}

// Set the manual wiener process mode
template <class C>
void Euler<C>::setManualWienerMode( const bool s ) { manualWiener=s; }
template <class C>
void Euler<C>::setWienerIncrements( const typename C::array a ) { dw=a; }

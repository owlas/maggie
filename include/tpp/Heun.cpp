// Heun.cpp
// Implementation for the Heun integration scheme
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//


// Constructor
template <class EQ>
Heun<EQ>::Heun( const EQ &le, const typename EQ::array& init_state,
               const double time, const double dt, mt19937& rng )
    : Integrator<EQ>( le, init_state, time)
    , h( dt )
    , dim( EQ::dim )
    , wDim( EQ::wdim )
    , dist(0,1)
    , gen( rng, dist )
{
    // empty
}

template <class C>
void Heun<C>::step()
{
    // Generate 1D array of Wiener increments'
    if( !manualWiener )
        for( auto &i : dw )
            i = sqrt(h) * gen();

    // PREDICTION
    this->getLE().computeDrift( tmp1, this->state, this->getTime() );
    this->getLE().computeDiffusion( tmp2, this->state, this->getTime() );
    for( unsigned int i=0; i!=dim; i++ )
    {
        xPred[i] = this->state[i] + tmp1[i]*h;
        for ( unsigned int j=0; j!=dim; j++ )
            xPred[i] += tmp2[i][j]*dw[j];
    }

    this->getLE().computeDrift( tmp1Up, this->xPred, this->getTime()+h );
    this->getLE().computeDiffusion( tmp2Up, this->xPred, this->getTime()+h );

    for( unsigned int i=0; i<dim; i++ )
    {
        xPred[i] = this->state[i] + 0.5*h*( tmp1Up[i] + tmp1[i] );
        for( unsigned int j=0; j<wDim; j++ )
            xPred[i] += 0.5*dw[j]*( tmp2Up[i][j] + tmp2[i][j] );
    }

    this->setState( xPred );
    this->setTime( this->getTime()+h );
}

// Set the manual wiener process mode
template <class C>
void Heun<C>::setManualWienerMode( const bool s ) { manualWiener=s; }
template <class C>
void Heun<C>::setWienerIncrements( const typename C::array& a ) { dw=a; }

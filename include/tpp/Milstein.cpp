// Milstein.cpp
// Implementation of the Milstein Taylor expansion numerical
// solver for stochastic differential equations. Implemented here for the
// stochastic Landau-Lifshitz-Gilbert equation in its reduced form
//
// Oliver W. Laslett (2015)
// O.Laslett@soton.ac.uk
//

// Constructor
template <class C>
Milstein<C>::Milstein( const C &le, const typename C::array& init_state,
                       const double time, const double dt,
                       mt19937& rng_1, mt19937& rng_2 )
    : Integrator<C>( le, init_state, time )
    , h( dt )
    , dim( C::dim )
    , wDim( C::wdim )
    , dist_1( 0,sqrt( h ) )
    , dist_2( 0,1 )
    , gen_1( rng_1, dist_1 )
    , gen_2( rng_2, dist_2 )
    , p( 200 )
    , eta( boost::extents[wDim][p] )
    , zeta( boost::extents[wDim][p] )
{
    // empty
}

template <class C>
void Milstein<C>::step()
{
    // Compute the drift and diffusion at the current time step
    this->getLE().computeDrift( tmp1, this->state, this->getTime() );
    this->getLE().computeDiffusion( tmp2, this->state, this->getTime() );
    this->getLE().computeDiffusionDerivatives( tmp3, this->state, this->getTime() );

    // Compute a Wiener process increment (unless manual mode)
    if( !manualWiener )
        for( auto &i : dw )
            i = gen_1();

    // Scaled Wiener process for computation of multiple integrals below
    for( unsigned int i=0; i<wDim; ++i )
        dw2[i] = dw[i]/sqrt( h );

    // rh0 for the calculation below
    double rho=0.0;
    for( unsigned int r=0; r!=p; ++r )
        rho += 1/pow( r+1, 2 );
    rho *= 1/( 2*pow( M_PI, 2 ) );
    rho += 1.0/12.;

    // Independent Gaussian RVs for integrals
    for( unsigned int i=0; i!=dim; i++ )
        *( mu.data()+i ) = gen_2();
    for( unsigned int i=0; i!=eta.num_elements(); i++ )
        *( eta.data()+i ) = gen_2();
    for( unsigned int i=0; i!=zeta.num_elements(); i++ )
        *( zeta.data()+i ) = gen_2();

    // Compute the multiple Stratonovich integrals
    for( unsigned int j1=0; j1!= wDim; ++j1 )
        for( unsigned int j2=0; j2!=wDim; ++j2 )
        {
            if( j1==j2 )
                tmp22[j1][j2] = 0.5*pow( dw[j1],2 );
            else
            {
                tmp22[j1][j2] = 0;
                for( bidx r=0; r!=p; ++r )
                    tmp22[j1][j2] += ( zeta[j1][r]*( sqrt( 2 )*dw2[j2]
                                                     + eta[j2][r] )
                                       - zeta[j2][r]*( sqrt( 2 )*dw2[j1]
                                                       + eta[j1][r] ) ) / (r+1);
                tmp22[j1][j2] *= h/( 2*M_PI );
                tmp22[j1][j2] += h*( 0.5*dw2[j1]*dw2[j2] + sqrt( rho )
                                     *( mu[j1]*dw2[j2] - mu[j2]*dw2[j1] ) );
            }
        }

    // Compute the next state for each vector
    for( unsigned int k=0; k!=dim; ++k )
    {
        next_state[k] = this->state[k] + tmp1[k]*h;

        for( unsigned int j=0; j!=wDim; j++ )
            next_state[k] += tmp2[k][j]*dw[j];

        for( unsigned int j1=0; j1!=wDim; j1++ )
            for( unsigned int j2=0; j2!=wDim; j2++ )
                for( unsigned int j3=0; j3!=dim; j3++ )
                    next_state[k] += tmp2[j3][j1]*tmp3[k][j2][j3]*tmp22[j1][j2];
    }

    this->setState( next_state );
    this->setTime( this->getTime()+h );
}

// Getters and setters for the truncation number
template <class C>
void Milstein<C>::setP( const int pset ) { p=pset; }
template <class C>
int Milstein<C>::getP() const { return p; }

// Set the manual wiener process mode
template <class C>
void Milstein<C>::setManualWienerMode( const bool s ) { manualWiener=s; }
template <class C>
void Milstein<C>::setWienerIncrements( const typename C::array a ) { dw=a; }

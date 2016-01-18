// Heun.cpp
// Implementation for the Heun integration scheme
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#include<Heun.hpp>


// Constructor
template <typename T>
Heun<T>::Heun( const SDE<T> &le, const array<T>& init_state,
               const T time, const T dt, mt19937& rng )
    : Integrator<SDE<T>,T>( le, init_state, time)
    , h( dt )
    , dim( le.getDim() )
    , wDim( le.getWDim() )
    , dw( boost::extents[wDim])
    , xPred( boost::extents[dim])
    , tmp1( boost::extents[dim])
    , tmp1Up( boost::extents[dim])
    , tmp2( boost::extents[dim][wDim] )
    , tmp2Up( boost::extents[dim][wDim] )
    , dist(0,1)
    , gen( rng, dist )
{
    // empty
}

template <typename T>
void Heun<T>::step()
{
    // Generate 1D array of Wiener increments
    if( !manualWiener )
        for( auto &i : dw )
            i = sqrt(h) * gen();

    // PREDICTION
    this->getLE().computeDrift( tmp1, this->state, this->getTime() );
    this->getLE().computeDiffusion( tmp2, this->state, this->getTime() );
    for( bidx i=0; i!=dim; i++ )
    {
        xPred[i] = this->state[i] + tmp1[i]*h;
        for ( bidx j=0; j!=dim; j++ )
            xPred[i] += tmp2[i][j]*dw[j];
    }

    this->getLE().computeDrift( tmp1Up, this->xPred, this->getTime()+h );
    this->getLE().computeDiffusion( tmp2Up, this->xPred, this->getTime()+h );

    for( int i=0; i<dim; i++ )
    {
        xPred[i] = this->state[i] + 0.5*h*( tmp1Up[i] + tmp1[i] );
        for( int j=0; j<wDim; j++ )
            xPred[i] += 0.5*dw[j]*( tmp2Up[i][j] + tmp2[i][j] );
    }

    this->setState( xPred );
    this->setTime( this->getTime()+h );
}

// Set the manual wiener process mode
template <typename T>
void Heun<T>::setManualWienerMode( const bool s ) { manualWiener=s; }
template <typename T>
void Heun<T>::setWienerIncrements( const array<T> a ) { dw=a; }

// Explicit template declaration
template class Heun<float>;
template class Heun<double>;

// RK4.cpp
// Implementation of the RK4 deterministic integrator
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//

// Constructor
template <class C>
RK4<C>::RK4( const C &le, const typename C::array& init_state,
             const double init_time, const double dt )
    : Integrator<C>( le, init_state, init_time )
    , h( dt )
    , dim( C::dim )
{}


template <class C>
void RK4<C>::step()
{
    // Step 1
    this->getLE().computeDrift( k1, this->state, this->getTime() );
    for( unsigned int i=0; i<dim; i++ )
        tmp[i] = k1[i]*h/2.0 + this->state[i];

    // Step 2
    this->getLE().computeDrift( k2, tmp, this->getTime() + h/2.0);
    for( unsigned int i=0; i<dim; i++ )
        tmp[i] = k2[i]*h/2.0 + this->state[i];

    // Step 3
    this->getLE().computeDrift( k3, tmp, this->getTime() + h/2.0 );
    for( unsigned int i=0; i<dim; i++ )
        tmp[i] = k3[i]*h + this->state[i];

    // Step 4 and update state
    this->getLE().computeDrift( k4, tmp, this->getTime() + h);
    for( unsigned int i=0; i<dim; i++ )
        tmp[i] = this->state[i]
            + ( h/6 )*( k1[i] + 2*k2[i] + 2*k3[i] + k4[i] );
    this->setState( tmp );

    // Update time
    this->setTime( this->getTime() + h );
}

// explicit template instantiation
// template class RK4<float>;
// template class RK4<double>;

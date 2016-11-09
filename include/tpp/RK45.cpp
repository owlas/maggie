// RK45.cpp
// Implementation of the adaptive RK45 deterministic integrator
//
// Oliver W. Laslett 2016
// O.Laslett@soton.ac.uk
//

// Constructor
template <class C>
RK45<C>::RK45( const C &le, const typename C::array& init_state,
               const double init_time, const double _eps )
    : Integrator<C>( le, init_state, init_time )
    , h( 1e-3 )
    , dim( C::dim )
    , eps( _eps )
{}


template <class C>
void RK45<C>::step()
{
    bool step_success = false;
    double err;
    while( step_success == false )
    {
        // k1
        this->getLE().computeDrift( k1, this->state, this->getTime() );

        // k2
        for( unsigned int i=0; i<dim; i++ )
            tmp[i] = this->state[i] + h*c11*k1[i];
        this->getLE().computeDrift( k2, tmp, this->getTime() + h*hc1 );

        // k3
        for( unsigned int i=0; i<dim; i++ )
            tmp[i] = this->state[i] + h*( c21*k1[i] + c22*k2[i] );
        this->getLE().computeDrift( k3, tmp, this->getTime() + h*hc2 );

        // k4
        for( unsigned int i=0; i<dim; i++ )
            tmp[i] = this->state[i]
                + h*( c31*k1[i] + c32*k2[i] + c33*k3[i] );
        this->getLE().computeDrift( k4, tmp, this->getTime() + h*hc3 );

        // k5
        for( unsigned int i=0; i<dim; i++ )
            tmp[i] = this->state[i]
                + h*( c41*k1[i]
                      + c42*k2[i]
                      + c43*k3[i]
                      + c44*k4[i]
                    );
        this->getLE().computeDrift( k5, tmp, this->getTime() + h*hc4 );

        // k6
        for( unsigned int i=0; i<dim; i++ )
            tmp[i] = this->state[i]
                + h*( + c51*k1[i]
                      + c52*k2[i]
                      + c53*k3[i]
                      + c54*k4[i]
                      + c55*k5[i]
                    );
        this->getLE().computeDrift( k6, tmp, this->getTime() + h*hc5 );

        // Compute order 5 estimate
        for( unsigned int i=0; i<dim; i++ )
            tmp[i] = this->state[i]
                + h*( x11*k1[i]
                      + x13*k3[i]
                      + x14*k4[i]
                      + x16*k6[i]
                    );

        // Compute order 4 estimate
        for( unsigned int i=0; i<dim; i++ )
            tmp2[i] = this->state[i]
                + h*( x21*k1[i]
                      + x23*k3[i]
                      + x24*k4[i]
                      + x25*k5[i]
                      + x26*k6[i]
                    );

        // Compute the error between the two approximations
        // scale according to ( eps + eps*|state| )
        err=0;
        double mag=0;
        for( unsigned int i=0; i<dim; i++ )
            mag = this->state[i] * this->state[i];
        mag = pow( mag, 0.5 );
        for( unsigned int i=0; i<dim; i++ )
            err += pow( std::abs( tmp[i] - tmp2[i] ), 2 );
        err = pow( err, 0.5 );
        err /= ( dim*eps*( 1 + mag ) );

        // If relative error is below 1 then step was successful
        // otherwise reduce the step size (max 10x reduction)
        if( err < 1.0 )
            step_success = true;
        else
        {
            double hfactor = 0.84*pow( err, -0.2 );
            hfactor = std::abs( hfactor ) < 0.1 ? 0.1 : hfactor;
            this->setStepSize( h*hfactor );
        }
    }

    // Set the new state, time, step size (max 5x increase)
    double hfactor = err==0.0 ? 5.0 : 0.84*pow( err, -0.2 );
    hfactor = hfactor>5 ? 5.0 : hfactor;
    this->setTime( this->getTime() + h );
    this->setState( tmp2 );
    this->setStepSize( h*hfactor );
}

template <class C>
void RK45<C>::setStepSize( const double _h ) { h=_h; }

// explicit template instantiation
// template class RK4<float>;
// template class RK4<double>;

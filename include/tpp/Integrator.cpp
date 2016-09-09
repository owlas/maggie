// Integrator.tpp
// Implementation for the abstract Integrator class
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//

// Constructor
template<class EQ_C>
Integrator<EQ_C>::Integrator( const EQ_C& le, const typename EQ_C::array& init_state, const double time )
: state( init_state )
    , lEq( &le )
    , initial_state( init_state )
    , t( time )
    , initial_t( time )
{}

// Get state
template<class EQ_C>
const typename EQ_C::array& Integrator<EQ_C>::getState() const { return state; }
// Set State
template<class EQ_C>
void Integrator<EQ_C>::setState( const typename EQ_C::array& s )
{
    state=s;
}

// get time
template<class EQ_C>
double Integrator<EQ_C>::getTime() const { return t; }
// set time
template<class EQ_C>
void Integrator<EQ_C>::setTime( const double time ) { t = time; }

// get the pointer to the integrator
template<class EQ_C>
const EQ_C& Integrator<EQ_C>::getLE() const { return *lEq; }

// reset the integrator to the initial condition
template<class EQ_C>
void Integrator<EQ_C>::reset()
{
  state = initial_state;
  t = initial_t;
}

// reset the integrator with a new inital condition
template<class EQ_C>
void Integrator<EQ_C>::reset( const typename EQ_C::array new_init )
{
  initial_state = new_init;
  reset();
}

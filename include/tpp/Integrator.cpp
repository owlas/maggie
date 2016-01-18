// Integrator.tpp
// Implementation for the abstract Integrator class
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//

// Constructor
template<class C, typename T>
Integrator<C, T>::Integrator( const C& le, const array<T>& init_state,
			   const T time )
  : state( boost::extents[le.getDim()] )
  , lEq( &le )
  , initial_state( boost::extents[le.getDim()] )
{
  initial_t = time;
  setTime( time );
  initial_state = init_state;
  setState( init_state );
}

// Get state
template<class C, typename T>
array<T> Integrator<C, T>::getState() const { return state; }
// Set State
template<class C, typename T>
void Integrator<C, T>::setState( const array<T> s )
{
if( int( s.shape()[0] ) != lEq->getDim() )
    {
      throw invalid_argument( "Error: state dimension does not match"
			      " Langevin equation dimension");
    }
  else
    {
      state = s;
    }
}

// get time
template<class C, typename T>
T Integrator<C, T>::getTime() const { return t; }
// set time
template<class C, typename T>
void Integrator<C, T>::setTime( const T time ) { t = time; }

// get the pointer to the integrator
template<class C, typename T>
const C& Integrator<C, T>::getLE() const { return *lEq; }

// reset the integrator to the initial condition
template<class C, typename T>
void Integrator<C, T>::reset()
{
  state = initial_state;
  t = initial_t;
}

// reset the integrator with a new inital condition
template<class C, typename T>
void Integrator<C, T>::reset( const array<T> new_init )
{
  initial_state = new_init;
  reset();
}


// Explicitly instantiate template instances
// Only these type combinations are visible to user
template class Integrator<ODE<float>, float>;
template class Integrator<ODE<float>, double>;
template class Integrator<SDE, float>;
template class Integrator<SDE, double>;

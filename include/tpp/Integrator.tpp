// Integrator.tpp
// Implementation for the abstract Integrator class
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//

// Constructor
template<class T>
Integrator<T>::Integrator( const T& le, const array_f& init_state,
			   const float time )
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
template<class T>
array_f Integrator<T>::getState() const { return state; }
// Set State
template<class T>
void Integrator<T>::setState( const array_f s )
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
template<class T>
float Integrator<T>::getTime() const { return t; }
// set time
template<class T>
void Integrator<T>::setTime( const float time ) { t = time; }

// get the pointer to the integrator
template<class T>
const T& Integrator<T>::getLE() const { return *lEq; }

// reset the integrator to the initial condition
template<class T>
void Integrator<T>::reset()
{
  state = initial_state;
  t = initial_t;
}

// reset the integrator with a new inital condition
template<class T>
void Integrator<T>::reset( const array_f new_init )
{
  initial_state = new_init;
  reset();
}


// Explicitly instantiate template instances
template class Integrator<ODE>;
template class Integrator<SDE>;

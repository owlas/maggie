// Simulation.cpp
// 
// Implementation of the simulation objects
// 
// Oliver W. Laslett (2015)
// O.Laslett@soton.ac.uk
// 
#include<Simulation.hpp>

// ----- CONSTRUCTOR -----
Simulation::Simulation( const ParticleCluster g, matrix_f init_state,
			float h, unsigned int n, float temp,
			array_f field )
  : geom( g )
  , state( boost::extents[geom.getNParticles()][3] )
  , stability( geom.computeStability( temp ) )
{
  setTimeStep( h );
  setSimLength( n );
  setState( init_state );
  setTemp( temp );
  setField( field );
}



// ----- Setter methods -----
void Simulation::setSimLength( unsigned int n ) { N=n; }

void Simulation::setTimeStep( float h ) { dt=h; }

void Simulation::setField( array_f field ) { h = field; }

void Simulation::setTemp( float t )
{
  if( t < 0 )
    throw std::invalid_argument("Temperature T must be greater than"
				" zero.");
  T = t;
}

void Simulation::setState( matrix_f s )
{
  size_t ndims=2, dim1=geom.getNParticles(), dim2=3;
  if( ( s.dimensionality != ndims ) or
      ( s.shape()[0] != dim1 ) or
      ( s.shape()[1] != dim2 ) )
    throw std::invalid_argument( "State of the system must be (Nx2)"
				 " - where N is the number of "
				 "particles" );
  state = s;
}



// ----- Getter methods -----
float Simulation::getTimeStep() const { return dt; }

unsigned int Simulation::getSimLength() const { return N; }

matrix_f Simulation::getState() const { return state; }

float Simulation::getTemp() const { return T; }

array_f Simulation::getField() const { return h; }

ParticleCluster Simulation::getGeometry() const
{ 
  return geom;
}


// ----- COMPUTE ARRHENIUS RATES -----
// determines the energy barriers and then the transition rate matrix
// from the exponential Neel-Arrhenius law
matrix_f Simulation::arrheniusMatrix() const
{
  throw "NOT IMPLEMENTED";
}
matrix_f Simulation::computeEnergyBarriers() const
{
  throw "NOT IMPLMENTED";
}

// Run the simulation - currently an integration of the Langevin
// dynamics without field
int run()
{
  // First you need a random generator
  mt19937 rng( 301 );

  // Then a particle cluster
  array_f h( boost::extents[3] );
  std::fill( h.data(), h.data()+h.num_elements(), 0 );
  std::vector<Particle> plist;
  Particle p1( 1, 2, 3, 4, 5, h );
  Particle p2( 2, 3, 2, 5, 1, h );
  plist.push_back( p1 );
  plist.push_back( p2 );
  matrix_f locs( boost::extents[2][3] );
  std::fill( locs.data(), locs.data() + locs.num_elements(), 0 );
  ParticleCluster clust( plist, locs );

  // Then set an initial condition

  // Compute field

  // Integrate every particle at each timestep

  return 1; // everything was fine
}

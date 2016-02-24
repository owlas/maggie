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
                        float stepsize, unsigned int n, float temp,
                        array_f field )
  : geom( g )
  , state( boost::extents[geom.getNParticles()][3] )
  , stability( geom.computeStability( temp ) )
  , h( boost::extents[3] )
{
  setTimeStep( stepsize );
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
int Simulation::run()
{
    // Currently only runs a single particle
    if( geom.getNParticles() != 1 )
        throw "Simulation only runs for a single particle.";

    // Get the particle info
    Particle p = geom.getParticle( 0 );

    // The thermal field intensity for reduced simulation
    float therm_strength{std::sqrt( p.getAlpha() * Constants::KB * T
                                     / ( p.getK() * p.getV()
                                         * ( 1 + std::pow( p.getAlpha(),2 ) ) ) ) };


    // Compute the reduced time for the simulation
    float Hk{ 2 * p.getK() / ( Constants::MU0 * p.getMs() ) };
    float t_factor{ Constants::GYROMAG * Constants::MU0 * Hk
            / ( 1+pow( p.getAlpha(), 2 ) ) };
    float dtau = dt * t_factor;


    // Set up an LLG equation for the particle
    // finite temperature with the reduced field
    StocLLG<float> llg( therm_strength, p.getAlpha(),
                        h[0]/Hk, h[1]/Hk, h[2]/Hk );


    // Get the initial condition of the particle
    array_f init( boost::extents[3] );
    for( array_f::index i=0; i!=3; ++i )
        init[i] = state[0][i];


    // Initialise the integrator and its RNG
    mt19937 rng( 301 );
    auto inte = Heun<float>( llg, init, 0.0, dtau, rng );


    // store the solutions here
    matrix_f sols( boost::extents[N][3] );


    // Step the integrator and store the result at each step
    for( array_f::index n=0; n!=N; ++n )
    {
        inte.step();
        // store the solutions
        for( array_f::index i=0; i!=3; ++i )
            sols[n][i] = inte.getState()[i];
    }


    // Write the results to the hardrive
    boostToFile( sols, "llgsim.out" );

  return 1; // everything was fine
}

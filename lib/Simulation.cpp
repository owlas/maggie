// Simulation.cpp
//
// Implementation of the simulation objects
//
// Oliver W. Laslett (2015)
// O.Laslett@soton.ac.uk
//
#include<Simulation.hpp>

// ----- CONSTRUCTOR -----
Simulation::Simulation( const ParticleCluster g, const ad_vec init_state,
                        const double stepsize, const unsigned int n,
                        const double temp, const array_d field )
  : geom( g )
  , h( boost::extents[3] )
  , stability( geom.computeStability( temp ) )
{
  setTimeStep( stepsize );
  setSimLength( n );
  setState( init_state );
  setTemp( temp );
  setField( field );
}



// ----- Setter methods -----
void Simulation::setSimLength( const unsigned int n ) { N=n; }

void Simulation::setTimeStep( const double h ) { dt=h; }

void Simulation::setField( const array_d field ) { h = field; }

void Simulation::setTemp( const double t )
{
  if( t < 0 )
    throw std::invalid_argument("Temperature T must be greater than"
                                " zero.");
  T = t;
}

void Simulation::setState( const ad_vec s )
{
  size_t dim1=geom.getNParticles(), dim2=3;
  if( ( s.size() != dim1 ) or
      ( s[0].shape()[0] != dim2 ) )
    throw std::invalid_argument( "State of the system must be (Nx2)"
                                 " - where N is the number of "
                                 "particles" );
  state = s;
}



// ----- Getter methods -----
double Simulation::getTimeStep() const { return dt; }

unsigned int Simulation::getSimLength() const { return N; }

const ad_vec& Simulation::getState() const { return state; }

double Simulation::getTemp() const { return T; }

array_d Simulation::getField() const { return h; }

ParticleCluster Simulation::getGeometry() const
{
  return geom;
}


// ----- COMPUTE ARRHENIUS RATES -----
// determines the energy barriers and then the transition rate matrix_d
// from the exponential Neel-Arrhenius law
matrix_d Simulation::arrheniusMatrix_D() const
{
  throw "NOT IMPLEMENTED";
}
matrix_d Simulation::computeEnergyBarriers() const
{
  throw "NOT IMPLMENTED";
}

// Run the simulation - currently an integration of the Langevin
// dynamics without field
int Simulation::runFull()
{
    // Currently only runs a single particle
    if( geom.getNParticles() != 1 )
        throw "Simulation only runs for a single particle.";

    // Get the particle info
    Particle p = geom.getParticle( 0 );

    // The thermal field intensity for reduced simulation
    double therm_strength{std::sqrt( p.getAlpha() * Constants::KB * T
                                     / ( p.getK() * p.getV()
                                         * ( 1 + std::pow( p.getAlpha(),2 ) ) ) ) };


    // Compute the reduced time for the simulation
    double Hk{ 2 * p.getK() / ( Constants::MU0 * p.getMs() ) };
    double t_factor{ Constants::GYROMAG * Constants::MU0 * Hk
            / ( 1+pow( p.getAlpha(), 2 ) ) };
    double dtau = dt * t_factor;


    // compute the reduced external field
    array_d happ( extents[3] );
    for( bidx i=0; i!=3; ++i )
        happ[i] = h[i]/Hk;

    // Set up an LLG equation for the particle
    // finite temperature with the reduced field
    StocLLG<double> llg( therm_strength, p.getAlpha(),
                        happ[0], happ[1], happ[2] );


    // Get the initial condition of the particle
    array_d init( boost::extents[3] );
    for( array_d::index i=0; i!=3; ++i )
        init[i] = state[0][i];


    // Initialise the integrator and its RNG
    mt19937 rng( 301 );
    auto inte = Heun<double>( llg, init, 0.0, dtau, rng );


    // store the solutions here
    matrix_d sols( boost::extents[N][3] );

    // reference to the current state and the field
    const array_d& currentState = inte.getState();
    array_d& currentField = llg.getReducedFieldRef();


    // Step the integrator and store the result at each step
    for( array_d::index n=0; n!=N; ++n )
    {
        // first compute the field due to the anisotropy
        // and put that in the field vector
        p.computeAnisotropyField( currentField, currentState );

        // then add the reduced external field
        for( bidx j=0; j!=3; ++j )
            currentField[j] += happ[j];

        // step the integrator
        inte.step();

        // store the solutions
        for( array_d::index i=0; i!=3; ++i )
            sols[n][i] = currentState[i];
    }


    // Write the results to the hardrive
    boostToFile( sols, "llg.mag" );

  return 1; // everything was fine
}

// Run the simulation
// this runs the system N times and stores the final state
// these states are then saved to the hard drive
int Simulation::runEnsemble( unsigned int Nruns )
{
    // Currently only runs a single particle
    if( geom.getNParticles() != 1 )
        throw "Simulation only runs for a single particle.";

    // Get the particle info
    Particle p = geom.getParticle( 0 );

    // The thermal field intensity for reduced simulation
    double therm_strength{std::sqrt( p.getAlpha() * Constants::KB * T
                                     / ( p.getK() * p.getV()
                                         * ( 1 + std::pow( p.getAlpha(),2 ) ) ) ) };


    // Compute the reduced time for the simulation
    double Hk{ 2 * p.getK() / ( Constants::MU0 * p.getMs() ) };
    double t_factor{ Constants::GYROMAG * Constants::MU0 * Hk
            / ( 1+pow( p.getAlpha(), 2 ) ) };
    double dtau = dt * t_factor;


    // compute the reduced external field
    array_d happ( extents[3] );
    for( bidx i=0; i!=3; ++i )
        happ[i] = h[i]/Hk;

    // Set up an LLG equation for the particle
    // finite temperature with the reduced field
    StocLLG<double> llg( therm_strength, p.getAlpha(),
                        happ[0], happ[1], happ[2] );


    // Get the initial condition of the particle
    array_d init( boost::extents[3] );
    for( array_d::index i=0; i!=3; ++i )
        init[i] = state[0][i];


    // Initialise the integrator and its RNG
    mt19937 rng( 301 );
    auto inte = Heun<double>( llg, init, 0.0, dtau, rng );


    // store the solutions here
    matrix_d sols( boost::extents[Nruns][3] );

    // reference to the current state and the field
    const array_d& currentState = inte.getState();
    array_d& currentField = llg.getReducedFieldRef();


    // Run the integrator a number of times
    for( bidx m=0; m!=Nruns; ++m )
    {
        // Reset the itegrator
        inte.reset();

        // Step the integrator and store the result at each step
        for( array_d::index n=0; n!=N; ++n )
        {
            // first compute the field due to the anisotropy
            // and put that in the field vector
            p.computeAnisotropyField( currentField, currentState );

            // then add the reduced external field
            for( bidx j=0; j!=3; ++j )
                currentField[j] += happ[j];

            // step the integrator
            inte.step();
        } // for each time step

        // Store the final state of the system
        for( bidx i=0; i!=3; ++i )
            sols[m][i] = currentState[i];
    }

    // Write the results to the hardrive
    boostToFile( sols, "ensemble.mag" );

  return 1; // everything was fine
}

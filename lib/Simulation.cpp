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
  equilibrium_rng.seed(1001);

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
    mt19937 rng2( 301 );
    auto inte = Heun<double>( llg, init, 0.0, dtau, rng2 );


    // store the solutions here
    matrix_d sols( boost::extents[N][3] );

    // reference to the current state and the field
    const array_d& currentState = inte.getState();
    array_d& currentField = llg.getReducedFieldRef();

    // temporary state vector
    array_d state_temp( boost::extents[3] );


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

        // renormalise the solution
        double norm = sqrt( currentState[0]*currentState[0]
                            + currentState[1]*currentState[1]
                            + currentState[2]*currentState[2] );
        for( bidx k=0; k!=3; ++k )
            state_temp[k] = currentState[k] / norm;
        inte.setState( state_temp );

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

// Computes the first passage time for a single particle in the absence of an
// external field
int Simulation::runFPT( const int N_ensemble, const bool alignup )
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

    // each simulation in the ensemble has its own RNG
    // generate seeds here
    mt19937 seed_rng( 8888 );
    std::vector<int> seeds;
    seeds.reserve( N_ensemble );
    boost::uniform_int<> int_dist;
    for( auto &i : seeds )
        i = int_dist( seed_rng );

    // store the fpt for each run
    array_d fpt( extents[N_ensemble] );

    // run for every simulation in the ensemble
    #pragma omp parallel for
    for( bidx i=0; i<N_ensemble; ++i )
    {
        // Set up an LLG equation for the particle
        // finite temperature with the reduced field
        StocLLG<double> llg( therm_strength, p.getAlpha(),
                             happ[0], happ[1], happ[2] );


        // The initial condition of the particle is determined by the flag
        // if align up is true, then particle is initially such that mz=1
        array_d init( boost::extents[3] );
        if( alignup )
        {
            init[0] = 0; init[1] = 0; init[2] = 1;
        }    // if it is false then the initial condition is drawn from the
        // equilibrium initial condition
        else
            init = equilibriumState();


        // Initialise the integrator and its RNG
        mt19937 heun_rng( seeds[i] );
        auto inte = Heun<double>( llg, init, 0.0, dtau, heun_rng );


        // reference to the current state and the field
        const array_d& currentState = inte.getState();
        array_d& currentField = llg.getReducedFieldRef();

        // Step the integrator until switch occurs
        // then store the time
        while( 1 )
        {
            // first compute the field due to the anisotropy
            // and put that in the field vector
            p.computeAnisotropyField( currentField, currentState );

            // then add the reduced external field
            for( bidx j=0; j!=3; ++j )
                currentField[j] += happ[j];

            // step the integrator
            inte.step();

            // check the end condition
            if( currentState[2] < 0 )
                break;
        }

        // store the time
        fpt[i] = inte.getTime();

    }

    // Write the results to the hardrive
    boostToFile( fpt, "llg.fpt" );

  return 1; // everything was fine
}


// returns a random state from the equilibrium distribution
array_d Simulation::equilibriumState()
{

    // Currently only runs a single particle
    if( geom.getNParticles() != 1 )
        throw "Simulation only runs for a single particle.";

    // Get the particle info
    Particle p = geom.getParticle( 0 );

    // define the pdf
    auto pdf = [&p, this]( double theta )
        {
            return sin( theta) * exp( -p.getK()*p.getV()*pow( sin( theta ), 2 )
                                      / ( Constants::KB * T ) );
        };

    // find the max
    array_d angles( extents[10000] );
    for( bidx i=0; i!=10000; ++i )
        angles[i] = i * M_PI_2 / 10000;
    double max = 0.0;
    for( auto &t : angles )
        if( pdf( t ) > max)
            max = pdf( t );

    // use uniform distribution as the candidate density
    // with M 10% higher that the max of target dist
    static boost::random::uniform_real_distribution<> gen_candidate( 0,M_PI );
    static boost::random::uniform_real_distribution<> gen_check( 0,1 );
    double M = 1.1 * max;


    // Rejection rate algorithm
    double candidate = 0;
    while( 1 )
    {
        candidate = gen_candidate( equilibrium_rng );
        double acceptance = pdf( candidate ) / M;
        double check = gen_check( equilibrium_rng );
        if( acceptance > check )
            break;
    }

    array_d init_state( extents[3] );
    init_state[0] = 0;
    init_state[1] = 0;
    init_state[2] = cos( candidate );

    return init_state;
}

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
                        const temperature temp, const field field,
                        const unsigned long int seed )
    : geom( g )
    , Nparticles( geom.getNParticles() )
    , sigmas( geom.computeStability( temp ) )
    , integrators_seed( seed )
{
    setTimeStep( stepsize );
    setSimLength( n );
    setState( init_state );
    setTemp( temp );
    setField( field );
    equilibrium_rng.seed(1001);

    // Assert that all particles have the same Ms value
    double ms( geom.getParticle( 0 ).getMs() );
    for( unsigned int i=0; i!=Nparticles; ++i )
        if( geom.getParticle( i ).getMs() != ms )
            throw "Simulation requires that all particles in the cluster have "
                "the same Ms value";


    // Set up an LLG equation for each particle in the ensemble
    // geom -- applied field reduced
    // per particle:
    //   thermal strength
    //   alpha
    std::vector<StocLLG> llgs;
    for( unsigned int i=0; i!=Nparticles; ++i )
    {
        // Get the particle
        auto p_ref = geom.getParticle( i );

        // The thermal field intensity for reduced simulation
        double therm_strength {
            std::sqrt( p_ref.getAlpha() * Constants::KB * T
                       / ( geom.getAverageAnisConstant() * p_ref.getV()
                           * ( 1 + std::pow( p_ref.getAlpha(),2 ) ) ) )
                };

        // Initialise the LLG equation
        maggie::field heff{ 0, 0, 0 };
        llgs.push_back(
            StocLLG( therm_strength, p_ref.getAlpha(), heff )
            );
    }


    // Compute the reduced field constant
    double Hk{
        2 * geom.getAverageAnisConstant()
            / ( Constants::MU0 * ms ) };


    // compute the reduced external field
    maggie::field happ{ h[0]/Hk, h[1]/Hk, h[2]/Hk };


    // Compute the reduced time for the simulation
    // TODO: check reduction factor
    // TODO: what is happening with this alpha reference?
    double t_factor{ Constants::GYROMAG * Constants::MU0 * Hk
            / ( 1+pow( geom.getParticle( 0 ).getAlpha(), 2 ) ) };
    reduced_time_factor = t_factor;
    double dtau = dt * t_factor;


    // init magnetic state controller
    // requires:
    //  geometry
    //  external field
    //  llg equation list
    //  integrators list
    simulationController = MagneticStateController<Euler<StocLLG>>::unique_ptr(
        new MagneticStateController<Euler<StocLLG>>(
            geom, happ, llgs, state, dtau, integrators_seed
            )
        );
}


// ----- Setter methods -----
void Simulation::setSimLength( const unsigned int n ) { Nsteps=n; }

void Simulation::setTimeStep( const double h ) { dt=h; }

void Simulation::setField( const field _h ) { h = _h; }

void Simulation::setTemp( const temperature _T )
{
  if( _T < 0 )
    throw std::invalid_argument("Temperature T must be greater than"
                                " zero.");
  T = _T;
}

void Simulation::setState( const std::vector<moment> s )
{
  if( s.size() != Nparticles )
    throw std::invalid_argument( "Vector of moments must be N."
                                 " - where N is the number of "
                                 "particles" );
  state = s;
}



// ----- Getter methods -----
double Simulation::getTimeStep() const { return dt; }

unsigned int Simulation::getSimLength() const { return Nsteps; }

const std::vector<moment>& Simulation::getState() const { return state; }

temperature Simulation::getTemp() const { return T; }

field Simulation::getField() const { return h; }

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

int Simulation::init()
{


    return 1;
}

// Run the simulation - currently an integration of the Langevin
// dynamics without field
int Simulation::runFull()
{

    // Initialise the simulation
    simulationController->reset();

    // store the solutions here
    std::vector<matrix_d> sols;
    for( unsigned int i=0; i!=Nparticles; ++i )
    {
        matrix_d sol_mat( boost::extents[Nsteps][3] );
        sols.push_back( sol_mat );
    }

    // Step the integrator and store the result at each step
    for( unsigned int n=0; n!=Nsteps; ++n )
    {
        // Step the system
        simulationController->step();

        // store the solutions
        std::vector<maggie::moment> state = simulationController->getState();
        for( unsigned int p=0; p!=Nparticles; ++p )
        {
            sols[p][n][0] = state[p][0];
            sols[p][n][1] = state[p][1];
            sols[p][n][2] = state[p][2];
        }

    }


    // Write the results to the hardrive
    for( unsigned int i=0; i!=Nparticles; ++i )
    {
        std::ostringstream fname;
        fname << "llg" << i << ".mag";
        boostToFile( sols[i], fname.str() );
    }

    // Compute the azimuth angles (cos m_z)
    for( unsigned int i=0; i!=Nparticles; ++i )
    {
        array_d azimuth( boost::extents[Nsteps] );
        for( unsigned int step=0; step!=Nsteps; ++step )
            azimuth[step] = std::acos( sols[i][step][2] );

        std::ostringstream fname;
        fname << "llg" << i << ".azimuth";
        boostToFile( azimuth, fname.str() );
    }

  return 1; // everything was fine
}

// Run the simulation
// this runs the system N times and stores the final state
// these states are then saved to the hard drive
int Simulation::runEnsemble( unsigned int Nruns )
{
    // store the solutions here
    matrix_d sols( boost::extents[Nruns][3] );

    // initialise the simulation
    simulationController->reset();

    // Run the integrator a number of times
    for( bidx m=0; m!=Nruns; ++m )
    {

        // Step the integrator and store the result at each step
        for( array_d::index n=0; n!=Nsteps; ++n )
        {
            // step the system
            simulationController->step();
        } // for each time step

        // Store the final state of the first particle
        auto currentstate = simulationController->getState();
        for( unsigned int i=0; i!=3; ++i )
            sols[m][i] = currentstate[0][i];

        // reset the integrators
        simulationController->reset();
    }

    // Write the results to the hardrive
    boostToFile( sols, "ensemble.mag" );

  return 1; // everything was fine
}

// Stores the magnetisation of the particle cluster
// for each in the ensemble
int Simulation::runFullEnsemble( unsigned int Nruns )
{
    // store the solutions here
    matrix_d sols( boost::extents[Nruns][Nsteps] );

    // initialise the simulation
    simulationController->reset();

    // Run the integrator a number of times
    for( bidx m=0; m<Nruns; ++m )
    {

        // Step the integrator and store the result at each step
        for( array_d::index n=0; n!=Nsteps; ++n )
        {
            // step the system
            simulationController->step();

            // Store the current state
            auto currentstate = simulationController->getState();

            // compute the total magnetisation
            sols[m][n] = 0;
            for( unsigned int p=0; p!=Nparticles; ++p )
                sols[m][n] += currentstate[p][2] / Nparticles;

        } // for each time step

        // reset the integrators
        simulationController->reset();
    }

    // Write the results to the hardrive
    std::ostringstream fname;
    fname << "fullEnsemble-process" << omp_get_thread_num() << ".mag";
    boostToFile( sols, fname.str() );

  return 1; // everything was fine
}


// Computes the first passage time for a single particle in the absence of an
// external field
int Simulation::runFPT( const int N_ensemble, const bool alignup )
{
    // Currently only runs a single particle
    if( geom.getNParticles() != 1 )
        throw "Simulation only runs for a single particle.";

    simulationController->reset();

    // store the fpt for each run
    array_d fpt( extents[N_ensemble] );

    // run for every simulation in the ensemble
    //#pragma omp parallel for schedule(dynamic, 1)
    for( bidx i=0; i<N_ensemble; ++i )
    {

        // The initial condition of the particle is determined by the flag
        // if align up is true, then particle is initially such that mz=1
        maggie::moment init;
        if( alignup )
        {
            init = maggie::moment{ 0, 0, 1 };
        }    // if it is false then the initial condition is drawn from the
        // equilibrium initial condition
        else
            init = equilibriumState()[0];

        std::vector<maggie::moment> init_states{ init };
        simulationController->reset( init_states );

        // Step the integrator until switch occurs
        // then store the time
        while( 1 )
        {
            // step the simulation
            simulationController->step();

            // check the end condition
            auto currentState = simulationController->getState()[0];
            if( currentState[2] < -0.5 )
                break;
        }

        // store the time
        fpt[i] = simulationController->getTime() / reduced_time_factor;
        cout << i << endl;

    }

    // Write the results to the hardrive
    auto stability_ratios = geom.computeStability( T );
    auto sr = stability_ratios[0];

    std::ostringstream fname;
    fname << "llg" << sr << ".fpt";
    boostToFile( fpt, fname.str() );

    return 1; // everything was fine
}


// returns a random state from the equilibrium distribution
std::vector<maggie::moment> Simulation::equilibriumState()
{

    std::vector<maggie::moment> random_state;

    // Get the particle info and compute the stability
    for( unsigned int p=0; p!=Nparticles; ++p )
    {
        // Define the pdf
        auto particle = geom.getParticle( p );
        double stability_ratio =
            particle.getK() * particle.getV() / ( Constants::KB * T );

        auto pdf = [stability_ratio]( double theta )
            {
                return sin( theta ) * exp( -stability_ratio * pow( sin( theta ), 2 ) );
            };

        // find the max
        array_d angles( extents[1000] );
        for( bidx i=0; i!=1000; ++i )
            angles[i] = i * M_PI_4 / 1000;
        double max = 0.0;
        for( auto &t : angles )
            if( pdf( t ) > max)
                max = pdf( t );

        // use uniform distribution as the candidate density
        // with M 10% higher that the max of target dist
        static boost::random::uniform_real_distribution<> gen_candidate( 0,M_PI_4 );
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

        random_state.push_back( maggie::moment{ 0, 0, cos( candidate ) } );
    }

    return random_state;
}

// Computes the time that switches occur for a single particle
int Simulation::runResidence( const unsigned int N_switches )
{
    // Currently only runs a single particle
    if( geom.getNParticles() != 1 )
        throw "Simulation only runs for a single particle.";

    // initialise
    simulationController->reset();

    // store the time of each switch
    array_d switch_times( extents[N_switches] );

    // Number of switches
    unsigned int n = 0;

    // Check if the initial state is up or down
    maggie::moment currentState = simulationController->getState()[0];
    int up;
    if ( currentState[2] > 0 )
        up = 1;
    else
        up = -1;

    // Step the integrator until switch occurs
    // then store the time
    while( n!=N_switches )
    {
        // step the simulation
        simulationController->step();

        // check the end condition
        currentState = simulationController->getState()[0];
        if( currentState[2]*up < 0 )
        {
            up *= -1;
            switch_times[n] = simulationController->getTime();
            n++;
            cout << n << endl;
        }
    }

    // Write the results to the hardrive
    boostToFile( switch_times, "llg.switches" );

    return 1; // everything was fine
}

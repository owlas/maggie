// llgsim.cpp - an example of using the simulation object to run a
//            - simulation of the llg dynamics of a single
//            - non interacting particle.
//
// Oliver W. Laslett (2016)
// O.Laslett@soton.ac.uk
//
#include <maggie.hpp>
#include <vector>
using plist=std::vector<Particle>;
#include <boost/multi_array.hpp>
using array_f=boost::multi_array<float, 1>;
using matrix_f=boost::multi_array<float, 2>;
using boost::extents;
using bidx=boost::multi_array_types::index;

int main()
{
    // The first thing to do is to set up the particle cluster
    // in this case we just have one particle, so lets define that:
    float gyromagnetic_constant{ Constants::GYROMAG };
    float damping{ 0.1 };
    float Ms{ 1400e3 };
    float diameter{ 5e-9 };
    float anisotropy_const{ 1 };
    array_f anisotropy_axis( extents[3] );
    anisotropy_axis[0] = 0.0;
    anisotropy_axis[1] = 0.0;
    anisotropy_axis[2] = 1.0;
    auto p = Particle( gyromagnetic_constant, damping, Ms, diameter,
                       anisotropy_const, anisotropy_axis );


    // The particle cluster then needs a list of particles and their
    // respective locations. In this case its easy.
    plist ps;
    ps.push_back( p );

    matrix_f locs( extents[1][3] );
    for( auto &i : locs[0] )
        i = 0.0;

    auto cluster = ParticleCluster( ps, locs );


    // Now we have defined the particle parameters and geometry we can
    // set up the environment

    // first the initial state of the magnetisation
    matrix_f init( extents[1][3] );
    init[0][0] = 1.0;
    init[0][1] = 0.0;
    init[0][2] = 0.0;

    // now the time step and simulation length
    float dt{ 1e-14 };
    int time_steps{ 100000 };

    // Finally we choose the temperature and external field
    float temp{ 30 };
    array_f field( extents[3] );
    field[0] = 0.0; field[1] = 0.0; field[2] = 100e3;

    auto mysim = Simulation( cluster, init, dt, time_steps, temp, field );


    // The run command runs an LLG simulation of the particle for the
    // specified number of time steps and saves the magnetisation data
    // to the hard disk
    mysim.run();


    // and that is that! open up pandas and take a look
    return 1;
}

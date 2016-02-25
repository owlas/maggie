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
using array_d=boost::multi_array<double, 1>;
using matrix_d=boost::multi_array<double, 2>;
using boost::extents;
using bidx=boost::multi_array_types::index;

#include <inih/INIReader.h>
#include <iostream>
using std::cout; using std::endl;
int main()
{
    // Read the ini file to get the particle properties
    // by default, read llgsim.ini
    INIReader reader( "llgsim.ini");
    if(reader.ParseError() < 0 )
    {
        cout << "Cannot load llgsim.ini" << endl;
        return 1;
    }

    // The first thing to do is to set up the particle cluster
    // in this case we just have one particle, so lets define that:
    double gyromagnetic_constant{ Constants::GYROMAG };
    double damping{ reader.GetReal( "particle", "alpha", -1 ) }; // 0.01
    double Ms{ reader.GetReal( "particle", "ms", -1 ) }; // 1400e3
    double diameter{ reader.GetReal( "particle", "diameter", -1 ) }; // 5e-9
    double anisotropy_const{ reader.GetReal( "particle", "anisotropy", -1 ) };
    array_d anisotropy_axis( extents[3] );
    anisotropy_axis[0] = reader.GetReal( "particle", "anis_x", -1 );
    anisotropy_axis[1] = reader.GetReal( "particle", "anis_y", -1 );
    anisotropy_axis[2] = reader.GetReal( "particle", "anis_z", -1 );
    auto p = Particle( gyromagnetic_constant, damping, Ms, diameter,
                       anisotropy_const, anisotropy_axis );


    // The particle cluster then needs a list of particles and their
    // respective locations. In this case its easy.
    plist ps;
    ps.push_back( p );

    matrix_d locs( extents[1][3] );
    for( auto &i : locs[0] )
        i = 0.0;

    auto cluster = ParticleCluster( ps, locs );


    // Now we have defined the particle parameters and geometry we can
    // set up the environment

    // first the initial state of the magnetisation
    array_d init( extents[3] );
    init[0] = 1.0;
    init[1] = 0.0;
    init[2] = 0.0;
    std::vector<array_d> states = { init };

    // now the time step and simulation length
    double dt{ 1e-14 };
    int time_steps{ 100000 };

    // Finally we choose the temperature and external field
    double temp{ 30 };
    array_d field( extents[3] );
    field[0] = 0.0; field[1] = 0.0; field[2] = 100e3;

    auto mysim = Simulation( cluster, states, dt, time_steps, temp, field );


    // The run command runs an LLG simulation of the particle for the
    // specified number of time steps and saves the magnetisation data
    // to the hard disk
    mysim.run();


    // and that is that! open up pandas and take a look
    return 1;
}

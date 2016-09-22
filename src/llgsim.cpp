// llgsim.cpp - an example of using the simulation object to run a
//            - simulation of the llg dynamics of a single
//            - non interacting particle.
//
// Oliver W. Laslett (2016)
// O.Laslett@soton.ac.uk
//
#include <maggie.hpp>
#include <vector>

#include <boost/multi_array.hpp>
using array_d=boost::multi_array<double, 1>;
using matrix_d=boost::multi_array<double, 2>;
using boost::extents;
using bidx=boost::multi_array_types::index;

#include "../include/types.hpp"

#include "../include/inih/INIReader.h"
#include <iostream>
using std::cout; using std::endl;
using std::flush;
int main()
{
    // Read the ini file to get the particle properties
    // by default, read llgsim.ini
    INIReader reader( "llgsim.ini");

    // A particle cluster needs a list of particles and locations
    // Initialise those vectors here.
    Particle::vector ps;
    std::vector<maggie::position> locs;
    std::vector<maggie::moment> states;

    if(reader.ParseError() < 0 )
    {
        cout << "Cannot load llgsim.ini" << endl;
        return 1;
    }

    // Global parameters
    double temp, dt;
    int time_steps;
    maggie::field field;

    // For each section of the ini file specifies a single particle
    for( auto &section_name : reader.GetSections() )
    {
      if( section_name == "global")
        {
          temp =  reader.GetReal( section_name, "temperature", -1 );
          dt = reader.GetReal( section_name, "time_step", -1 );
          time_steps = reader.GetInteger( section_name, "number_steps", -1 );
          field = maggie::field{
              reader.GetReal( section_name, "Happ_x", -1 ),
              reader.GetReal( section_name, "Happ_y", -1 ),
              reader.GetReal( section_name, "Happ_z", -1 )
            };
        }
      else
        {
          maggie::magnetogyric gyromagnetic_constant{ Constants::GYROMAG };
          maggie::damping damping{ reader.GetReal( section_name, "alpha", -1 ) };
          maggie::magnetisation Ms{ reader.GetReal( section_name, "ms", -1 ) };
          maggie::diameter diameter{ reader.GetReal( section_name, "diameter", -1 ) };
          maggie::anisotropy anisotropy_const{ reader.GetReal( section_name, "anisotropy", -1 ) };
          maggie::axis anisotropy{
            reader.GetReal( section_name, "anis_x", -1 ),
              reader.GetReal( section_name, "anis_y", -1 ),
              reader.GetReal( section_name, "anis_z", -1 )
              };
          auto p = Particle( gyromagnetic_constant, damping, Ms, diameter,
                             anisotropy_const, anisotropy );
          maggie::position loc{
            reader.GetReal( section_name, "loc_x", -1 ),
              reader.GetReal( section_name, "loc_y", -1 ),
              reader.GetReal( section_name, "loc_z", -1 )
              };
          maggie::moment state{
            reader.GetReal( section_name, "m_x", -1 ),
              reader.GetReal( section_name, "m_y", -1 ),
              reader.GetReal( section_name, "m_z", -1 )
              };
          states.push_back( state );
          ps.push_back( p );
          locs.push_back( loc );
        }
    }


    // Create the cluster
    auto cluster = ParticleCluster( ps, locs );

    // print stability
    for( unsigned int p=0; p!=cluster.getNParticles(); ++p )
        cout << "stability ratio: " << ps[p].getK() * ps[p].getV() / ( Constants::KB * temp ) << endl;

    auto mysim = Simulation( cluster, states, dt, time_steps, temp, field );

    // The run command runs an LLG simulation of the particle for the
    // specified number of time steps and save the time domain signal
    // of the magnetisation to the hard disk
    cout << "Computing example trajectory... " << flush;
    mysim.runFull();
    cout << "done" << endl;


    // Run an ensemble of the system and save the final state of the each
    // member of the ensemble
    cout << "Computing equilibrium distribution... " << flush;
    mysim.runEnsemble( 10000 );
    cout << "done" << endl;

    // and that is that! open up pandas and take a look
    return 1;
}

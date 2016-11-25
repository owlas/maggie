#ifndef MAGTRAJ_H
#define MAGTRAJ_H

#include "types.hpp"
#include <stdio.h>
#include <boost/multi_array.hpp>
#include <string>
#include <memory>

class MagnetisationTrajectory
{
public:
    MagnetisationTrajectory( const int n_particles, const int n_time_steps );
    int save_to_file( std::string fname, bool compress=false ) const;
    boost::multi_array<double,3>& get_trajectories_ref();
    boost::multi_array<double,1>& get_time_ref();
    boost::multi_array<double,3>::reference operator[]( const int index );

    // Useful aliases
    using u_ptr = std::unique_ptr<MagnetisationTrajectory>;

private:
    const size_t n_particles;
    const size_t n_time_steps;
    const size_t write_batch_size=200;
    size_t number_of_writes;
    size_t remainder_batch;
    boost::multi_array<double,3> trajectories;
    boost::multi_array<double,1> time;
};
#endif

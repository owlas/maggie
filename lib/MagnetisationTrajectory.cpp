#include "../include/MagnetisationTrajectory.hpp"
#include <zlib.h>

MagnetisationTrajectory::MagnetisationTrajectory(
    const int _n_particles, const int _n_time_steps )
    : n_particles( _n_particles )
    , n_time_steps( _n_time_steps )
    , trajectories( boost::extents[n_time_steps][n_particles][3] )
    , time( boost::extents[n_time_steps] )
{
    number_of_writes = n_time_steps / write_batch_size;
    remainder_batch = n_time_steps % write_batch_size;
}

boost::multi_array<double,3> &MagnetisationTrajectory::get_trajectories_ref()
{
    return trajectories;
}

boost::multi_array<double,1> &MagnetisationTrajectory::get_time_ref()
{
    return time;
}

int MagnetisationTrajectory::save_to_file( std::string fname, bool compress ) const
{
    if( compress )
    {
        gzFile out;
        void *ptr;
        out = gzopen( fname.c_str(), "wb" );
        if( out==NULL ) return 1;
        size_t offset;
        for( unsigned int i=0; i!=number_of_writes; ++i )
        {
            offset = 3*n_particles*i*write_batch_size;
            ptr = (void*) ( trajectories.origin() + offset );
            gzwrite( out, ptr, sizeof(double) *write_batch_size );
        }
        offset = 3*n_particles*number_of_writes*write_batch_size;
        ptr = (void*) ( trajectories.origin() + offset );
        gzwrite( out, ptr, sizeof(double) * remainder_batch * 3 * n_particles );
        gzclose( out );
    }
    else
    {
        FILE* out;
        out = fopen( fname.c_str(), "wb" );
        if( out==NULL ) return 1;
        size_t offset;
        for( unsigned int i=0; i!=number_of_writes; ++i )
        {
            offset = 3*n_particles*i*write_batch_size;
            fwrite( trajectories.origin() + offset, sizeof(double), write_batch_size * 3 * n_particles, out );
        }
        offset = 3*n_particles*number_of_writes*write_batch_size;
        fwrite( trajectories.origin() + offset, sizeof(double), remainder_batch * 3 * n_particles, out );
        fclose( out );
    }

    return 0;
}

boost::multi_array<double,3>::reference MagnetisationTrajectory::operator[]( const int index )
{
    return trajectories[index];
}

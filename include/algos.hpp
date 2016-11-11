// Algorithms for use in maggie package
//
// Oliver W. Laslett (2016)
//
#ifndef ALGOS_H
#define ALGOS_H

#include <cstdlib>
namespace algos
{
    template<class Container>
    void interpolate_to_regular_line( Container const& values,
                                      Container const& grid_points,
                                      Container& regular_grid_values )
    {
        // calculate regular grid spacing
        double spacing = grid_points.back() / ( regular_grid_values.size() - 1 );

        // add start and end conditions
        regular_grid_values.front() = values.front();
        regular_grid_values.back() = values.back();

        // track how many points have been written
        unsigned int points_written=1;

        // map each point in the irregular grid onto the regular grid
        for( unsigned int i=1; i!=values.size(); ++i )
            while( points_written*spacing <= grid_points[i] )
            {
                regular_grid_values[points_written] = values[i-1]
                    + ( points_written*spacing - grid_points[i-1] )
                    / ( grid_points[i] - grid_points[i-1] )
                    * ( values[i] - values[i-1] );
                points_written++;
            }
    }
}

#endif

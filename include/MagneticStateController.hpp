#ifndef MAGCONTH
#define MAGCONTH


#include "types.hpp"
using namespace maggie;

#include<memory>
#include "ParticleCluster.hpp"
#include "StocLLG.hpp"
#include "Heun.hpp"
#include "Euler.hpp"
#include <iostream>

#include<boost/random.hpp>
using boost::mt19937;

#include <boost/multi_array.hpp>
using matrix_d = boost::multi_array<double, 2>;
using array3_d = boost::multi_array<double, 3>;

template<class INTE>
class MagneticStateController
{
public:
    // constructor
    MagneticStateController(
        const ParticleCluster geom,
        const maggie::field appliedField,
        const StocLLG::vector,
        const std::vector<maggie::moment> init_states,
        const double time_step,
        const unsigned int seed
        );

    // useful aliases
    using unique_ptr = std::unique_ptr<MagneticStateController>;

    // compute the effective field
    void updateEffectiveField();

    // Renormalise the field lengths
    void renormaliseStates();

    // Compute one step
    void step();

    // Reset state to initial
    void reset();

    // reset to state
    void reset( std::vector<maggie::moment> );

    // getters
    double getDDPrefactor() const;
    matrix_d getCubeDisplacements() const;
    std::vector<moment> getState() const;
    double getTime() const;

private:
    ParticleCluster geom;
    maggie::field happ;
    StocLLG::vector llgs;
    int N;
    array3_d disp;
    matrix_d cube_disps;
    std::vector<maggie::volume> v;
    double ddNorm;
    mt19937 integrator_rng;
    typename INTE::vector integrators;
    std::vector<maggie::field*> heff;
};

#include "tpp/MagneticStateController.cpp"
#endif

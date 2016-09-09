// StocLLG.hpp
// Header for the stochastic LLG Langevin equation
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef STOCLLG_H
#define STOCLLG_H

#include "types.hpp"
#include "SDE.hpp"
#include <vector>

class StocLLG : public SDE<3,3>
{
 public:
    // Constructor when the thermal intensity (s) and damping (a) is
    // calculated in advance
    StocLLG( const maggie::stability s, const maggie::damping a,
             const maggie::field h );

    // useful aliases
    using vector = std::vector<StocLLG>;

    // returns the drift vector
    virtual void computeDrift( StocLLG::array&, const StocLLG::array&, const double ) const;
    // returns the diffusion matrix
    virtual void computeDiffusion( StocLLG::matrix&, const StocLLG::array&,
                                   const double ) const;
    // returns a vector of derivative sums for Taylor
    virtual void computeDiffusionDerivatives( StocLLG::array3 &out, const StocLLG::array &in,
                                              const double ) const;

    void setReducedHeff( maggie::field );
    maggie::field getReducedHeff() const;

    void setSigma( const maggie::stability );
    maggie::stability getSigma() const;

    void setAlpha( const maggie::damping );
    maggie::damping getAlpha() const;

    // get a modifiable reference to the field
    maggie::field& getReducedFieldRef();

 private:
    maggie::field h;
    maggie::stability sigma;
    maggie::damping alpha;
};

// This is the Ito version of the Stochastic LLG
class StocLLGIto : public StocLLG
{
public:
    StocLLGIto( const maggie::stability s, const maggie::damping a,
                const maggie::field h );

    // returns the drift vector
    void computeDrift( StocLLGIto::array&, const StocLLGIto::array&, const double ) const;

};

#endif

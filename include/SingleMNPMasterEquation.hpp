// SingleMNPMasterEquation.hpp - Class to represent the master
// equation for a single uniaxial magnetic nanoparticle. Derived from
// LangevunEquation base.
//
// Oliver Laslett (2015)
// O.Laslett@soton.ac.uk
//
#include<LangevinEquation.hpp>
#include<Constants.hpp>
#include<Field.hpp>
#include<stdexcept>
#include<functional>

class SingleMNPMasterEquation : public LangevinEquation
{
public:
  // constructor
  SingleMNPMasterEquation( float temperature,
                           float anisotropy_constant,
                           float paticle_radius,
                           std::function<float( float )> field_func,
                           float field_angle );

  // compute the difference in energy between the two states of the
  // system
  float ediff( float t );  

  // Compute the drift of the master equation
  virtual void computeDrift( array_f& out, array_f& in, float t );

private:
  float k;
  float T;
  float r;
  float v;
  float psi;
  std::function<float(float)> field;
  const float TAU0=10e-10;
};

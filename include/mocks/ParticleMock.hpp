#include <gmock.hpp>
#include "../Particle.hpp"

class MockParticle : public Particle
{
    MOCK_CONST_METHOD0( getK, double() );
    MOCK_CONST_METHOD0( getV, double() );
}

#ifndef TYPES_H
#define TYPES_H

#include <array>

namespace maggie {
    using magnetogyric = double;
    using damping = double;
    using diameter = double;
    using anisotropy = double;
    using magnetisation = double;
    using volume = double;
    using temperature = double;
    using stability = double;
    using axis = std::array<double,3>;
    using field = std::array<double,3>;
    using moment = std::array<double,3>;
    using position = std::array<double,3>;
}

#endif

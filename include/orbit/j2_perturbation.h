#pragma once

#include "orbit/body.h"

namespace orbit {

class J2Perturbation {
public:
  // earth и sat заданы в одной инерциальной системе координат, в AU.
  [[nodiscard]] Vec3 accelerationEarthJ2(const Body &earth,
                                         const Body &sat) const;
};

} // namespace orbit

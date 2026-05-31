#include "orbit/j2_perturbation.h"

#include "orbit/constants.h"

#include <cmath>
#include <stdexcept>

namespace orbit {

Vec3 J2Perturbation::accelerationEarthJ2(const Body &earth,
                                         const Body &sat) const {
  const double x = sat.r[0] - earth.r[0];
  const double y = sat.r[1] - earth.r[1];
  const double z = sat.r[2] - earth.r[2];

  const double r2 = x * x + y * y + z * z;

  if (r2 == 0.0) {
    throw std::runtime_error("acceleration_J2_earth: zero distance");
  }

  const double r = std::sqrt(r2);
  const double r5 = r2 * r2 * r;

  const double z2_r2 = z * z / r2;

  const double k = 1.5 * constants::EarthJ2 * earth.GM *
                   constants::EarthRadiusAu * constants::EarthRadiusAu / r5;

  return {k * x * (5.0 * z2_r2 - 1.0), k * y * (5.0 * z2_r2 - 1.0),
          k * z * (5.0 * z2_r2 - 3.0)};
}

} // namespace orbit

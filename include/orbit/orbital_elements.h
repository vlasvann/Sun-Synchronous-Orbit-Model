#pragma once

#include "orbit/body.h"

namespace orbit {

struct OrbitalElements {
  double semiMajorAxisAu{};           // a, AU
  double eccentricity{};              // e
  double inclinationRad{};            // i, radians
  double longitudeAscendingNodeRad{}; // Omega, radians
  double argumentPericenterRad{};     // omega, radians
  double trueAnomalyRad{};            // nu, radians
};

class OrbitalElementsCalculator {
public:
  // earth и satellite должны быть заданы в одной системе координат.
  // Для элементов орбиты используются относительные r = r_sat - r_earth
  // и v = v_sat - v_earth. mu = earth.GM.
  [[nodiscard]] static OrbitalElements fromState(const Body &earth,
                                                 const Body &satellite);
};

}
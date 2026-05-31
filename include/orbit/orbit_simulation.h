#pragma once

#include "orbit/body.h"
#include "orbit/rk8_integrator.h"

#include <iosfwd>
#include <vector>

namespace orbit {

class OrbitSimulation {
public:
  explicit OrbitSimulation(const Rk8Integrator &integrator);

  void run(std::vector<Body> &bodies, std::ostream &output) const;

  [[nodiscard]] static double distanceKm(const Body &a, const Body &b);
  [[nodiscard]] static double relativeSpeedKmS(const Body &a, const Body &b);

private:
  const Rk8Integrator &integrator_;
};

} // namespace orbit

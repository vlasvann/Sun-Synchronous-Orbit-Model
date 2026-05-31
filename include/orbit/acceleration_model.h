#pragma once

#include "orbit/body.h"
#include "orbit/j2_perturbation.h"

#include <vector>

namespace orbit {

class AccelerationModel {
public:
  // порядок тел сохраняется; acc[i] соответствует bodies[i].
  [[nodiscard]] std::vector<Vec3>
  computeNewtonAndRelativistic(const std::vector<Body> &bodies) const;

private:
  J2Perturbation j2_;
};

} // namespace orbit

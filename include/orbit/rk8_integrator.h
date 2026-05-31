#pragma once

#include "orbit/acceleration_model.h"
#include "orbit/body.h"

#include <vector>

namespace orbit {

class Rk8Integrator {
public:
  explicit Rk8Integrator(const AccelerationModel &accelerationModel);

  // возвращает новое состояние через шаг h.
  [[nodiscard]] std::vector<Body> step(const std::vector<Body> &bodies,
                                       double h) const;

private:
  const AccelerationModel &accelerationModel_;
};

} // namespace orbit

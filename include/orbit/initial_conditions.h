#pragma once

#include "orbit/body.h"

#include <vector>

namespace orbit {

class InitialConditionsFactory {
public:
  // Возвращает систему, используемую в расчётах: Earth + SSO.
  [[nodiscard]] static std::vector<Body> createSunSynchronousBodies();
};

} // namespace orbit

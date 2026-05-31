#pragma once

#include <array>
#include <string>

namespace orbit {

using Vec3 = std::array<double, 3>;

struct Body {
  std::string name;
  double GM{};
  Vec3 r{};
  Vec3 v{};
};

} // namespace orbit

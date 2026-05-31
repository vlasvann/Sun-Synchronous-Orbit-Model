#pragma once

#include "orbit/body.h"

#include <cmath>

namespace orbit::math {

inline double dot3(const Vec3 &a, const Vec3 &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline Vec3 cross3(const Vec3 &a, const Vec3 &b) {
  return {a[1] * b[2] - a[2] * b[1],
          a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

inline double norm3(const Vec3 &a) { return std::sqrt(dot3(a, a)); }

inline Vec3 add3(const Vec3 &a, const Vec3 &b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

inline Vec3 sub3(const Vec3 &a, const Vec3 &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

inline Vec3 mul3(const Vec3 &a, double factor) {
  return {factor * a[0], factor * a[1], factor * a[2]};
}

inline Vec3 div3(const Vec3 &a, double factor) {
  return {a[0] / factor, a[1] / factor, a[2] / factor};
}

} 
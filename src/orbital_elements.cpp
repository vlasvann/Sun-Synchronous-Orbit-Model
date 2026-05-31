#include "orbit/orbital_elements.h"

#include "orbit/constants.h"
#include "orbit/vector_math.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace orbit {
namespace {

constexpr double Eps = 1.0e-14;

[[nodiscard]] double clampUnit(double x) {
  return std::max(-1.0, std::min(1.0, x));
}

[[nodiscard]] double wrapToTwoPi(double angle) {
  const double twoPi = 2.0 * constants::Pi;
  angle = std::fmod(angle, twoPi);
  if (angle < 0.0) {
    angle += twoPi;
  }
  return angle;
}

} // namespace

OrbitalElements OrbitalElementsCalculator::fromState(const Body &earth,
                                                     const Body &satellite) {
  const Vec3 r = math::sub3(satellite.r, earth.r);
  const Vec3 v = math::sub3(satellite.v, earth.v);

  const double rNorm = math::norm3(r);
  const double vNorm = math::norm3(v);

  if (rNorm <= Eps) {
    throw std::runtime_error("OrbitalElementsCalculator: zero relative radius");
  }
  if (earth.GM <= 0.0) {
    throw std::runtime_error(
        "OrbitalElementsCalculator: central GM must be positive");
  }

  const double mu = earth.GM;

  // h = r x v - удельный орбитальный момент импульса.
  const Vec3 h = math::cross3(r, v);
  const double hNorm = math::norm3(h);
  if (hNorm <= Eps) {
    throw std::runtime_error(
        "OrbitalElementsCalculator: zero angular momentum");
  }

  // n = k x h - вектор линии узлов; k = ось Z базовой плоскости.
  const Vec3 kAxis = {0.0, 0.0, 1.0};
  const Vec3 n = math::cross3(kAxis, h);
  const double nNorm = math::norm3(n);

  // e_vec = (v x h) / mu - r / |r| - вектор эксцентриситета.
  const Vec3 eccentricityVector =
      math::sub3(math::div3(math::cross3(v, h), mu), math::div3(r, rNorm));
  const double e = math::norm3(eccentricityVector);

  const double specificEnergy = 0.5 * vNorm * vNorm - mu / rNorm;
  const double a = -mu / (2.0 * specificEnergy);

  OrbitalElements elements;
  elements.semiMajorAxisAu = a;
  elements.eccentricity = e;
  elements.inclinationRad = std::acos(clampUnit(h[2] / hNorm));

  // долгота восходящего узла
  if (nNorm > Eps) {
    elements.longitudeAscendingNodeRad = std::acos(clampUnit(n[0] / nNorm));
    if (n[1] < 0.0) {
      elements.longitudeAscendingNodeRad =
          2.0 * constants::Pi - elements.longitudeAscendingNodeRad;
    }
  } else {
    // Экваториальная орбита: линия узлов не определена.
    elements.longitudeAscendingNodeRad = 0.0;
  }

  if (nNorm > Eps && e > Eps) {
    elements.argumentPericenterRad =
        std::acos(clampUnit(math::dot3(n, eccentricityVector) / (nNorm * e)));
    if (eccentricityVector[2] < 0.0) {
      elements.argumentPericenterRad =
          2.0 * constants::Pi - elements.argumentPericenterRad;
    }
  } else {
    // Круговая или экваториальная орбита: аргумент перицентра вырожден.
    elements.argumentPericenterRad = 0.0;
  }

  if (e > Eps) {
    elements.trueAnomalyRad =
        std::acos(clampUnit(math::dot3(eccentricityVector, r) / (e * rNorm)));
    if (math::dot3(r, v) < 0.0) {
      elements.trueAnomalyRad = 2.0 * constants::Pi - elements.trueAnomalyRad;
    }
  } else {
    elements.trueAnomalyRad = 0.0;
  }

  elements.longitudeAscendingNodeRad =
      wrapToTwoPi(elements.longitudeAscendingNodeRad);
  elements.argumentPericenterRad = wrapToTwoPi(elements.argumentPericenterRad);
  elements.trueAnomalyRad = wrapToTwoPi(elements.trueAnomalyRad);

  return elements;
}

} // namespace orbit
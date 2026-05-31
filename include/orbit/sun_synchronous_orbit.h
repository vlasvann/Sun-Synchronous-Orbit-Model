#pragma once

#include "orbit/body.h"

namespace orbit {

struct KeplerianElementsInput {
  double semiMajorAxisKm{};           // a, km
  double eccentricity{};              // e
  double inclinationRad{};            // i, radians
  double longitudeAscendingNodeRad{}; // Omega, radians
  double argumentPericenterRad{};     // omega, radians
  double trueAnomalyRad{};            // nu, radians
};

class SunSynchronousOrbit {
public:
  // Теоретическая J2-прецессия линии узлов:
  // Omega_dot = -3/2 * n * J2 * (R_E / p)^2 * cos(i)
  // Возвращает rad/s.
  [[nodiscard]] static double nodalPrecessionRadS(double semiMajorAxisKm,
                                                  double eccentricity,
                                                  double inclinationRad);

  // f(i) = Omega_dot(i) - n_E.
  // Для SSO надо f(i) = 0.
  [[nodiscard]] static double sunSyncEquation(double semiMajorAxisKm,
                                              double eccentricity,
                                              double inclinationRad);

  // Метод бисекции по наклонению.
  // leftDeg и rightDeg задают отрезок поиска в градусах.
  [[nodiscard]] static double
  findInclinationBisectionRad(double semiMajorAxisKm, double eccentricity,
                              double leftDeg = 90.0, double rightDeg = 110.0,
                              double toleranceRad = 1.0e-13,
                              int maxIterations = 200);

  // Перевод кеплеровских элементов спутника относительно Земли
  // в барицентрические координаты тела.
  [[nodiscard]] static Body
  makeSatelliteFromKeplerian(const Body &earth,
                             const KeplerianElementsInput &elements,
                             const std::string &satelliteName = "SSO");
};

} // namespace orbit
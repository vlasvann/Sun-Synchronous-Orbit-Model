#include "orbit/sun_synchronous_orbit.h"

#include "orbit/constants.h"

#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <iostream>

namespace orbit {
namespace {

[[nodiscard]] double degToRad(double deg) { return deg * constants::DegToRad; }

[[nodiscard]] double sign(double x) {
  if (x > 0.0) {
    return 1.0;
  }
  if (x < 0.0) {
    return -1.0;
  }
  return 0.0;
}

} // namespace

double SunSynchronousOrbit::nodalPrecessionRadS(double semiMajorAxisKm,
                                                double eccentricity,
                                                double inclinationRad) {
  if (semiMajorAxisKm <= 0.0) {
    throw std::runtime_error(
        "nodalPrecessionRadS: semi-major axis must be positive");
  }

  if (eccentricity < 0.0 || eccentricity >= 1.0) {
    throw std::runtime_error(
        "nodalPrecessionRadS: eccentricity must be in [0, 1)");
  }

  const double a = semiMajorAxisKm;
  const double e = eccentricity;
  const double i = inclinationRad;

  const double p = a * (1.0 - e * e);

  const double n0 = std::sqrt(constants::EarthMuKm3S2 / (a * a * a));
  
  const double rOverP = constants::EarthRadiusKm / p;
  const double sinI = std::sin(i);

  const double nBar =
      n0 * (1.0 + 0.75 * constants::EarthJ2 * rOverP * rOverP *
                       std::sqrt(1.0 - e * e) *
                       (2.0 - 3.0 * sinI * sinI));

  return -1.5 * constants::EarthJ2 * rOverP * rOverP * nBar *
         std::cos(i);
}

double SunSynchronousOrbit::sunSyncEquation(double semiMajorAxisKm,
                                            double eccentricity,
                                            double inclinationRad) {
  return nodalPrecessionRadS(semiMajorAxisKm, eccentricity, inclinationRad) -
         constants::EarthMeanMotionAroundSunRadS;
}

double SunSynchronousOrbit::findInclinationBisectionRad(
    double semiMajorAxisKm, double eccentricity, double leftDeg,
    double rightDeg, double toleranceRad, int maxIterations) {
  double left = degToRad(leftDeg);
  double right = degToRad(rightDeg);

  double fLeft = sunSyncEquation(semiMajorAxisKm, eccentricity, left);
  double fRight = sunSyncEquation(semiMajorAxisKm, eccentricity, right);

  if (sign(fLeft) == 0.0) {
    return left;
  }

  if (sign(fRight) == 0.0) {
    return right;
  }

  if (sign(fLeft) == sign(fRight)) {
    throw std::runtime_error(
        "findInclinationBisectionRad: f(left) and f(right) have the same sign");
  }

  for (int iter = 0; iter < maxIterations; ++iter) {
    const double mid = 0.5 * (left + right);
    const double fMid = sunSyncEquation(semiMajorAxisKm, eccentricity, mid);
    std::cout << std::fixed << std::setprecision(12);

    std::cout << iter << " "
              << left * constants::RadToDeg << " "
              << right * constants::RadToDeg << " "
              << mid * constants::RadToDeg << " "
              << nodalPrecessionRadS(semiMajorAxisKm, eccentricity, mid) *
                    constants::RadToDeg * constants::SecondsPerDay << " "
              << fMid * constants::RadToDeg * constants::SecondsPerDay
              << "\n";

    if (std::abs(fMid) < 1.0e-16 || 0.5 * (right - left) < toleranceRad) {
      return mid;
    }

    if (sign(fLeft) != sign(fMid)) {
      right = mid;
      fRight = fMid;
      (void)fRight;
    } else {
      left = mid;
      fLeft = fMid;
    }
  }

  return 0.5 * (left + right);
}

Body SunSynchronousOrbit::makeSatelliteFromKeplerian(
    const Body &earth, const KeplerianElementsInput &elements,
    const std::string &satelliteName) {
  const double a = elements.semiMajorAxisKm;
  const double e = elements.eccentricity;
  const double i = elements.inclinationRad;
  const double Omega = elements.longitudeAscendingNodeRad;
  const double omega = elements.argumentPericenterRad;
  const double nu = elements.trueAnomalyRad;

  if (a <= 0.0) {
    throw std::runtime_error(
        "makeSatelliteFromKeplerian: semi-major axis must be positive");
  }

  if (e < 0.0 || e >= 1.0) {
    throw std::runtime_error(
        "makeSatelliteFromKeplerian: eccentricity must be in [0, 1)");
  }

  const double mu = constants::EarthMuKm3S2;
  const double p = a * (1.0 - e * e);

  const double rKm = p / (1.0 + e * std::cos(nu));

  // Координаты и скорость в перицентральной системе координат.
  const double xPf = rKm * std::cos(nu);
  const double yPf = rKm * std::sin(nu);

  const double velocityFactor = std::sqrt(mu / p);
  const double vxPf = -velocityFactor * std::sin(nu);
  const double vyPf = velocityFactor * (e + std::cos(nu));

  const double cosO = std::cos(Omega);
  const double sinO = std::sin(Omega);
  const double cosi = std::cos(i);
  const double sini = std::sin(i);
  const double cosw = std::cos(omega);
  const double sinw = std::sin(omega);

  // Матрица поворота R3(Omega) * R1(i) * R3(omega).
  const double r11 = cosO * cosw - sinO * sinw * cosi;
  const double r12 = -cosO * sinw - sinO * cosw * cosi;

  const double r21 = sinO * cosw + cosO * sinw * cosi;
  const double r22 = -sinO * sinw + cosO * cosw * cosi;

  const double r31 = sinw * sini;
  const double r32 = cosw * sini;

  const Vec3 rGeoKm = {
      r11 * xPf + r12 * yPf,
      r21 * xPf + r22 * yPf,
      r31 * xPf + r32 * yPf,
  };

  const Vec3 vGeoKmS = {
      r11 * vxPf + r12 * vyPf,
      r21 * vxPf + r22 * vyPf,
      r31 * vxPf + r32 * vyPf,
  };

  Body satellite;
  satellite.name = satelliteName;
  satellite.GM = 0.0;

  for (int k = 0; k < 3; ++k) {
    satellite.r[k] = earth.r[k] + rGeoKm[k] / constants::AuKm;
    satellite.v[k] =
        earth.v[k] + vGeoKmS[k] * constants::SecondsPerDay / constants::AuKm;
  }

  return satellite;
}

} // namespace orbit
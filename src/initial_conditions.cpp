#include "orbit/initial_conditions.h"
#include "orbit/constants.h"
#include "orbit/sun_synchronous_orbit.h"

#include <iostream>

namespace orbit {

std::vector<Body> InitialConditionsFactory::createSunSynchronousBodies() {
  Body earth;
  earth.name = "Earth";
  earth.GM = 8.887692445125e-10; // AU^3/day^2

  earth.r = {-7.814927467557258E-01, -5.904284172621328E-01,
             -2.558079827416129E-01};

  earth.v = {1.062424825330070E-02, -1.226562008833967E-02,
             -5.316635391581201E-03};

  // Типичный диапазон SSO: 600--800 km.
  const double altitudeKm = 800.0;
  const double semiMajorAxisKm = constants::EarthRadiusKm + altitudeKm;

  const double eccentricity = 0;

  const double inclinationRad =
      SunSynchronousOrbit::findInclinationBisectionRad(
          semiMajorAxisKm, eccentricity, 90.0, 110.0);

  const double nodalRateRadS = SunSynchronousOrbit::nodalPrecessionRadS(
      semiMajorAxisKm, eccentricity, inclinationRad);

  std::cout << "Sun-synchronous orbit parameters:\n";
  std::cout << "  altitude       = " << altitudeKm << " km\n";
  std::cout << "  a              = " << semiMajorAxisKm << " km\n";
  std::cout << "  e              = " << eccentricity << "\n";
  std::cout << "  i              = " << inclinationRad * constants::RadToDeg
            << " deg\n";
  std::cout << "  Omega_dot      = "
            << nodalRateRadS * constants::RadToDeg * constants::SecondsPerDay
            << " deg/day\n";
  std::cout << "  target n_E     = "
            << constants::EarthMeanMotionAroundSunDegDay << " deg/day\n";

  KeplerianElementsInput elements;
  elements.semiMajorAxisKm = semiMajorAxisKm;
  elements.eccentricity = eccentricity;
  elements.inclinationRad = inclinationRad;

  // Эти три угла выбираем произвольно.
  // Они задают ориентацию и начальную точку спутника на орбите.
  elements.longitudeAscendingNodeRad = 0.0 * constants::DegToRad;
  elements.argumentPericenterRad = 0.0 * constants::DegToRad;
  elements.trueAnomalyRad = 0.0 * constants::DegToRad;

  Body sso =
      SunSynchronousOrbit::makeSatelliteFromKeplerian(earth, elements, "SSO");

  return {earth, sso};
}

} // namespace orbit

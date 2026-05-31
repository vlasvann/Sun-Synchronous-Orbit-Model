#include "orbit/orbit_simulation.h"

#include "orbit/constants.h"
#include "orbit/orbital_elements.h"

#include <cmath>
#include <iomanip>
#include <ostream>

namespace orbit {

OrbitSimulation::OrbitSimulation(const Rk8Integrator &integrator)
    : integrator_(integrator) {}

double OrbitSimulation::distanceKm(const Body &a, const Body &b) {
  const double dx = b.r[0] - a.r[0];
  const double dy = b.r[1] - a.r[1];
  const double dz = b.r[2] - a.r[2];

  const double r_au = std::sqrt(dx * dx + dy * dy + dz * dz);
  return r_au * constants::AuKm;
}

double OrbitSimulation::relativeSpeedKmS(const Body &a, const Body &b) {
  const double dvx = b.v[0] - a.v[0];
  const double dvy = b.v[1] - a.v[1];
  const double dvz = b.v[2] - a.v[2];

  const double v_au_day = std::sqrt(dvx * dvx + dvy * dvy + dvz * dvz);
  return v_au_day * constants::AuKm / constants::SecondsPerDay;
}

void OrbitSimulation::run(std::vector<Body> &bodies,
                          std::ostream &output) const {
  output << std::setprecision(15);

  output << "# t_days jd_tdb distance_km relative_speed_km_s "
            "a_km eccentricity inclination_deg Omega_deg omega_deg "
            "true_anomaly_deg\n";

  for (int i = 0; i <= constants::SimulationSteps; ++i) {
    const double tDays = i * constants::TimeStepDays;
    const double jdTdb = constants::InitialJulianDateTdb + tDays;

    const double r = distanceKm(bodies[0], bodies[1]);
    const double v = relativeSpeedKmS(bodies[0], bodies[1]);

    const OrbitalElements elements =
        OrbitalElementsCalculator::fromState(bodies[0], bodies[1]);

    output << tDays << " " << jdTdb << " " << r << " " << v << " "
           << elements.semiMajorAxisAu * constants::AuKm << " "
           << elements.eccentricity << " "
           << elements.inclinationRad * constants::RadToDeg << " "
           << elements.longitudeAscendingNodeRad * constants::RadToDeg << " "
           << elements.argumentPericenterRad * constants::RadToDeg << " "
           << elements.trueAnomalyRad * constants::RadToDeg << "\n";

    bodies = integrator_.step(bodies, constants::TimeStepDays);
  }
}

} // namespace orbit
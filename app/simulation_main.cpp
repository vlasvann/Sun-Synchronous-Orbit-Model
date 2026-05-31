#include "orbit/acceleration_model.h"
#include "orbit/constants.h"
#include "orbit/initial_conditions.h"
#include "orbit/orbit_simulation.h"
#include "orbit/rk8_integrator.h"

#include <fstream>
#include <iostream>

int main() {
  auto bodies = orbit::InitialConditionsFactory::createSunSynchronousBodies();

  std::cout << "Initial Earth-satellite distance = "
            << orbit::OrbitSimulation::distanceKm(bodies[0], bodies[1])
            << " km\n";

  std::cout << "Initial Earth-SSO relative speed = "
            << orbit::OrbitSimulation::relativeSpeedKmS(bodies[0], bodies[1])
            << " km/s\n";

  std::ofstream file(orbit::constants::OutputFileName);

  if (!file) {
    std::cerr << "Cannot open orbit_output.txt\n";
    return 1;
  }

  const orbit::AccelerationModel accelerationModel;
  const orbit::Rk8Integrator integrator(accelerationModel);
  const orbit::OrbitSimulation simulation(integrator);

  simulation.run(bodies, file);

  std::cout << "Final Earth-SSO distance = "
            << orbit::OrbitSimulation::distanceKm(bodies[0], bodies[1])
            << " km\n";

  return 0;
}

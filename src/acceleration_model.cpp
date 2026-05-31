#include "orbit/acceleration_model.h"

#include "orbit/constants.h"
#include "orbit/vector_math.h"

#include <cmath>
#include <stdexcept>

namespace orbit {

std::vector<Vec3> AccelerationModel::computeNewtonAndRelativistic(
    const std::vector<Body> &bodies) const {
  const std::size_t n = bodies.size();
  if (n == 0) {
    throw std::runtime_error("compute_accelerations_newton: empty bodies list");
  }

  // acc0 содержит только ньютоновское ускорение.
  std::vector<Vec3> acc0(n, {0.0, 0.0, 0.0});
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      if (i == j) {
        continue;
      }

      const double dx = bodies[j].r[0] - bodies[i].r[0];
      const double dy = bodies[j].r[1] - bodies[i].r[1];
      const double dz = bodies[j].r[2] - bodies[i].r[2];
      const double r2 = dx * dx + dy * dy + dz * dz;
      const double r = std::sqrt(r2);
      if (r == 0.0) {
        throw std::runtime_error("zero distance");
      }

      const double inv_r3 = 1.0 / (r2 * r);
      acc0[i][0] += bodies[j].GM * dx * inv_r3;
      acc0[i][1] += bodies[j].GM * dy * inv_r3;
      acc0[i][2] += bodies[j].GM * dz * inv_r3;
    }
  }

  std::vector<Vec3> acc(n, {0.0, 0.0, 0.0});

  for (std::size_t A = 0; A < n; ++A) {
    acc[A] = acc0[A];

    const auto &vA = bodies[A].v;
    const double vA2 = math::dot3(vA, vA);

    for (std::size_t B = 0; B < n; ++B) {
      if (B == A) {
        continue;
      }

      const auto &vB = bodies[B].v;

      const double dx = bodies[B].r[0] - bodies[A].r[0];
      const double dy = bodies[B].r[1] - bodies[A].r[1];
      const double dz = bodies[B].r[2] - bodies[A].r[2];
      const double r2 = dx * dx + dy * dy + dz * dz;
      const double r = std::sqrt(r2);
      const double inv_r = 1.0 / r;
      const double inv_r2 = inv_r * inv_r;

      const Vec3 n_AB = {dx * inv_r, dy * inv_r, dz * inv_r};
      const Vec3 n_BA = {-dx * inv_r, -dy * inv_r, -dz * inv_r};

      const double vB2 = math::dot3(vB, vB);
      const double vA_vB = math::dot3(vA, vB);
      const double nAB_vB = math::dot3(n_AB, vB);

      double sumAC = 0.0;
      for (std::size_t C = 0; C < n; ++C) {
        if (C == A) {
          continue;
        }

        const double ex = bodies[C].r[0] - bodies[A].r[0];
        const double ey = bodies[C].r[1] - bodies[A].r[1];
        const double ez = bodies[C].r[2] - bodies[A].r[2];
        sumAC += bodies[C].GM / std::sqrt(ex * ex + ey * ey + ez * ez);
      }

      double sumBC = 0.0;
      for (std::size_t C = 0; C < n; ++C) {
        if (C == B) {
          continue;
        }

        const double ex = bodies[C].r[0] - bodies[B].r[0];
        const double ey = bodies[C].r[1] - bodies[B].r[1];
        const double ez = bodies[C].r[2] - bodies[B].r[2];
        sumBC += bodies[C].GM / std::sqrt(ex * ex + ey * ey + ez * ez);
      }

      const Vec3 rBA = {dx, dy, dz};
      const double rBA_aB = math::dot3(rBA, acc0[B]);

      const double S = vA2 + 2.0 * vB2 - 4.0 * vA_vB - 1.5 * nAB_vB * nAB_vB -
                       4.0 * sumAC - sumBC + 0.5 * rBA_aB;

      const double f2 =
          bodies[B].GM * inv_r2 / constants::SpeedOfLightAuPerDaySquared * S;
      acc[A][0] += f2 * n_BA[0];
      acc[A][1] += f2 * n_BA[1];
      acc[A][2] += f2 * n_BA[2];

      const Vec3 v4A3B = {4 * vA[0] - 3 * vB[0], 4 * vA[1] - 3 * vB[1],
                          4 * vA[2] - 3 * vB[2]};
      const double f3 = bodies[B].GM * inv_r2 /
                        constants::SpeedOfLightAuPerDaySquared *
                        math::dot3(n_AB, v4A3B);
      acc[A][0] += f3 * (vA[0] - vB[0]);
      acc[A][1] += f3 * (vA[1] - vB[1]);
      acc[A][2] += f3 * (vA[2] - vB[2]);

      const double f4 =
          3.5 * bodies[B].GM * inv_r / constants::SpeedOfLightAuPerDaySquared;
      acc[A][0] += f4 * acc0[B][0];
      acc[A][1] += f4 * acc0[B][1];
      acc[A][2] += f4 * acc0[B][2];
    }
  }

  // J2 добавляется для любого безмассового спутника вокруг Земли.
  // Предполагаем, что bodies[0] = Earth, bodies[1] = satellite.
  if (bodies.size() >= 2 && bodies[0].name == "Earth" && bodies[1].GM == 0.0) {
    const auto aJ2 = j2_.accelerationEarthJ2(bodies[0], bodies[1]);

    acc[1][0] += aJ2[0];
    acc[1][1] += aJ2[1];
    acc[1][2] += aJ2[2];
  }

  return acc;
}

} // namespace orbit

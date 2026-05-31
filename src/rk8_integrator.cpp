#include "orbit/rk8_integrator.h"

#include "orbit/constants.h"

#include <array>
#include <cstddef>

namespace orbit {

namespace {

constexpr std::size_t S = constants::Rk8StageCount;

// Коэффициенты RK8
constexpr double A[S][S] = {
    {0.0},

    {1.0 / 18.0},

    {1.0 / 48.0, 1.0 / 16.0},

    {1.0 / 32.0, 0.0, 3.0 / 32.0},

    {5.0 / 16.0, 0.0, -75.0 / 64.0, 75.0 / 64.0},

    {3.0 / 80.0, 0.0, 0.0, 3.0 / 16.0, 3.0 / 20.0},

    {29443841.0 / 614563906.0, 0.0, 0.0, 77736538.0 / 692538347.0,
     -28693883.0 / 1125000000.0, 23124283.0 / 1800000000.0},

    {16016141.0 / 946692911.0, 0.0, 0.0, 61564180.0 / 158732637.0,
     22789713.0 / 633445777.0, 545815736.0 / 2771057229.0,
     -180193667.0 / 1043307555.0},

    {39632708.0 / 573591083.0, 0.0, 0.0, -433636366.0 / 683701615.0,
     -421739975.0 / 2616292301.0, 100302831.0 / 723423059.0,
     790204164.0 / 839813087.0, 800635310.0 / 3783071287.0},

    {246121993.0 / 1340847787.0, 0.0, 0.0, -37695042795.0 / 15268766246.0,
     -309121744.0 / 1061227803.0, -12992083.0 / 490766935.0,
     6005943493.0 / 2108947869.0, 393006217.0 / 1396673457.0,
     123872331.0 / 1001029789.0},

    {-1028468189.0 / 846180014.0, 0.0, 0.0, 8478235783.0 / 508512852.0,
     1311729495.0 / 1432422823.0, -10304129995.0 / 1701304382.0,
     -48777925059.0 / 3047939560.0, 15336726248.0 / 1032824649.0,
     -45442868181.0 / 3398467696.0, 3065993473.0 / 597172653.0},

    {185892177.0 / 718116043.0, 0.0, 0.0, -3185094517.0 / 667107341.0,
     -477755414.0 / 1098053517.0, -703635378.0 / 230739211.0,
     5731566787.0 / 1027545527.0, 5232866602.0 / 850066563.0,
     -4093664535.0 / 808688257.0, 3962137247.0 / 1805957418.0,
     65686358.0 / 487910083.0},

    {403863854.0 / 491063109.0, 0.0, 0.0, -5068492393.0 / 434740067.0,
     -411421997.0 / 543043805.0, 652783627.0 / 914296604.0,
     11173962825.0 / 925320556.0, -13158990841.0 / 6184727034.0,
     3936647629.0 / 1978049680.0, -160528059.0 / 685178525.0,
     248638103.0 / 1413531060.0, 0.0}};

constexpr double B[S] = {14005451.0 / 335480064.0,
                         0.0,
                         0.0,
                         0.0,
                         0.0,
                         -59238493.0 / 1068277825.0,
                         181606767.0 / 758867731.0,
                         561292985.0 / 797845732.0,
                         -1041891430.0 / 1371343529.0,
                         760417239.0 / 1151165299.0,
                         118820643.0 / 751138087.0,
                         -528747749.0 / 2220607170.0,
                         1.0 / 4.0};

} // namespace

Rk8Integrator::Rk8Integrator(const AccelerationModel &accelerationModel)
    : accelerationModel_(accelerationModel) {}

std::vector<Body> Rk8Integrator::step(const std::vector<Body> &bodies,
                                      double h) const {
  const std::size_t n = bodies.size();

  std::array<std::vector<Vec3>, S> k_r;
  std::array<std::vector<Vec3>, S> k_v;

  for (std::size_t s = 0; s < S; ++s) {
    k_r[s].resize(n);
    k_v[s].resize(n);
  }

  for (std::size_t i = 0; i < n; ++i) {
    k_r[0][i] = bodies[i].v;
  }
  k_v[0] = accelerationModel_.computeNewtonAndRelativistic(bodies);

  for (std::size_t s = 1; s < S; ++s) {
    std::vector<Body> tmp = bodies;

    for (std::size_t i = 0; i < n; ++i) {
      for (int dim = 0; dim < constants::SpaceDimension; ++dim) {
        double sum_r = 0.0;
        double sum_v = 0.0;

        for (std::size_t j = 0; j < s; ++j) {
          sum_r += A[s][j] * k_r[j][i][dim];
          sum_v += A[s][j] * k_v[j][i][dim];
        }

        tmp[i].r[dim] += h * sum_r;
        tmp[i].v[dim] += h * sum_v;
      }
    }

    for (std::size_t i = 0; i < n; ++i) {
      k_r[s][i] = tmp[i].v;
    }

    k_v[s] = accelerationModel_.computeNewtonAndRelativistic(tmp);
  }

  std::vector<Body> out = bodies;

  for (std::size_t i = 0; i < n; ++i) {
    for (int dim = 0; dim < constants::SpaceDimension; ++dim) {
      double sum_r = 0.0;
      double sum_v = 0.0;

      for (std::size_t s = 0; s < S; ++s) {
        sum_r += B[s] * k_r[s][i][dim];
        sum_v += B[s] * k_v[s][i][dim];
      }

      out[i].r[dim] += h * sum_r;
      out[i].v[dim] += h * sum_v;
    }
  }

  return out;
}

} // namespace orbit

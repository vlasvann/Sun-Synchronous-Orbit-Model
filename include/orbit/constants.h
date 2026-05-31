#pragma once

#include <cstddef>

namespace orbit::constants {

inline constexpr double AuKm = 149597870.700;
inline constexpr double SecondsPerDay = 86400.0;

inline constexpr double Pi = 3.141592653589793238462643383279502884;
inline constexpr double RadToDeg = 180.0 / Pi;
inline constexpr double DegToRad = Pi / 180.0;

inline constexpr double SpeedOfLightAuPerDay = 173.144632674;
inline constexpr double SpeedOfLightAuPerDaySquared =
    SpeedOfLightAuPerDay * SpeedOfLightAuPerDay;

inline constexpr double EarthMuKm3S2 = 398600.435436;
inline constexpr double EarthRadiusKm = 6378.137;

inline constexpr double EarthRadiusAu = EarthRadiusKm / AuKm;
inline constexpr double EarthJ2 = 1.08262545e-3;

inline constexpr double EarthMeanMotionAroundSunRadS =
    0.9856473598947981 * DegToRad / SecondsPerDay;

inline constexpr double EarthMeanMotionAroundSunDegDay =
    EarthMeanMotionAroundSunRadS * RadToDeg * SecondsPerDay;

inline constexpr std::size_t Rk8StageCount = 13;
inline constexpr int SpaceDimension = 3;

inline constexpr double InitialJulianDateTdb = 2461160.500000000;

inline constexpr double TimeStepDays = 60.0 / SecondsPerDay;

inline constexpr int SimulationSteps = 90 * 24 * 60;

inline constexpr const char *OutputFileName = "../scripts/orbit_output.txt";

} // namespace orbit::constants
import numpy as np
import matplotlib.pyplot as plt

MU_EARTH_KM3_S2 = 398600.435436
R_EARTH_KM = 6378.137
J2_EARTH = 1.08262545e-3

SECONDS_PER_DAY = 86400.0
RAD_TO_DEG = 180.0 / np.pi
YEAR_DAYS = 365.0

SUN_RATE_DEG_DAY = 0.9856473598947981


def rad_s_to_deg_day(value_rad_s):
    """Перевод угловой скорости из рад/с в град/сутки."""
    return value_rad_s * RAD_TO_DEG * SECONDS_PER_DAY


def relative_error_percent(numerical, theoretical):
    """Относительная ошибка в процентах."""
    return abs(numerical - theoretical) / abs(theoretical) * 100.0


data = np.loadtxt("orbit_output.txt", comments="#")

t = data[:, 0]  # время, дни
distance_km = data[:, 2]  # расстояние от Земли до спутника, км
a_km = data[:, 4]  # большая полуось, км
eccentricity = data[:, 5]  # эксцентриситет
inclination_deg = data[:, 6]  # наклонение, градусы
Omega_deg = data[:, 7]  # долгота восходящего узла, градусы


mask_year = t <= YEAR_DAYS

t_year = t[mask_year]
distance_km_year = distance_km[mask_year]
a_km_year = a_km[mask_year]
eccentricity_year = eccentricity[mask_year]
inclination_deg_year = inclination_deg[mask_year]
Omega_deg_year = Omega_deg[mask_year]


# Убираем скачки угла через 360 градусов
Omega_unwrapped = np.rad2deg(np.unwrap(np.deg2rad(Omega_deg_year)))

# Линейная аппроксимация Omega(t)
Omega_fit = np.polyfit(t_year, Omega_unwrapped, 1)
Omega_dot_num = Omega_fit[0]
Omega_fit_values = np.polyval(Omega_fit, t_year)


# Идеальное изменение Omega для солнечно-синхронной орбиты
Omega_sun_like = Omega_unwrapped[0] + SUN_RATE_DEG_DAY * (t_year - t_year[0])

Omega_minus_sun_deg = Omega_unwrapped - Omega_sun_like

# 1 градус по долготе соответствует примерно 4 минутам LTAN
ltan_drift_min = 4.0 * Omega_minus_sun_deg


print()
print("Sun-synchronism check using LTAN drift, 1 year")
print(f"target Sun rate  = {SUN_RATE_DEG_DAY:.10f} deg/day")
print(f"Omega_dot_num    = {Omega_dot_num:.10f} deg/day")
print(f"LTAN drift min   = {np.min(ltan_drift_min):.6f} min")
print(f"LTAN drift max   = {np.max(ltan_drift_min):.6f} min")
print(f"LTAN drift range = {np.ptp(ltan_drift_min):.6f} min")


a_mean = np.mean(a_km_year)
e_mean = np.mean(eccentricity_year)
i_mean_deg = np.mean(inclination_deg_year)
i_mean_rad = np.deg2rad(i_mean_deg)

p_mean = a_mean * (1.0 - e_mean * e_mean)
n_mean = np.sqrt(MU_EARTH_KM3_S2 / (a_mean**3))


Omega_dot_theory_rad_s = (
    -1.5 * n_mean * J2_EARTH * (R_EARTH_KM / p_mean) ** 2 * np.cos(i_mean_rad)
)

Omega_dot_theory = rad_s_to_deg_day(Omega_dot_theory_rad_s)
Omega_error = relative_error_percent(Omega_dot_num, Omega_dot_theory)


print()
print("Mean orbital elements, 1 year")
print(f"a_mean           = {a_mean:.6f} km")
print(f"e_mean           = {e_mean:.10f}")
print(f"i_mean           = {i_mean_deg:.6f} deg")
print(f"p_mean           = {p_mean:.6f} km")
print(f"n_mean           = {n_mean:.10e} rad/s")

print()
print("Distance from Earth to satellite, 1 year")
print(f"distance_min     = {np.min(distance_km_year):.6f} km")
print(f"distance_max     = {np.max(distance_km_year):.6f} km")
print(f"distance_mean    = {np.mean(distance_km_year):.6f} km")

print()
print("Numerical nodal precession from linear fit")
print(f"Omega_dot_num    = {Omega_dot_num:.8f} deg/day")

print()
print("Theoretical J2 nodal precession")
print(f"Omega_dot_theory = {Omega_dot_theory:.8f} deg/day")

print()
print("Difference")
print(f"Omega error      = {Omega_error:.4f} %")


# дрейф LTAN за 1 год

plt.figure()
plt.plot(t_year, ltan_drift_min, linewidth=0.9, label="LTAN drift")
plt.axhline(0.0, linestyle="--", linewidth=1)

y_min = 0
y_max = 25
padding = 0.15 * (y_max - y_min)
plt.ylim(y_min - padding, y_max + padding)

plt.xlabel("t, дни")
plt.ylabel("дрейф LTAN, минуты")
plt.title("Проверка солнечной синхронизации: дрейф LTAN, 1 год")
plt.grid(True)
plt.legend()
plt.savefig("LTAN_drift_1_year.png", dpi=200)


# Прецессия Omega за 1 год

plt.figure()
plt.plot(t_year, Omega_unwrapped, label="Omega numerical")
plt.plot(t_year, Omega_fit_values, "--", label="linear fit")
plt.xlabel("t, дни")
plt.ylabel("Omega, градусы")
plt.title("Узловая прецессия: Omega(t), 1 год")
plt.grid(True)
plt.legend()
plt.savefig("Omega_precession_1_year.png", dpi=200)


# Расстояние Земля-спутник за первые 2 дня

mask_2_days = t_year <= 2.0

plt.figure()
plt.plot(
    t_year[mask_2_days],
    distance_km_year[mask_2_days],
    label="Расстояние между Землей и спутником",
)
plt.xlabel("t, дни")
plt.ylabel("расстояние, км")
plt.title("Расстояние от Земли до спутника, первые 2 дня")
plt.grid(True)
plt.legend()
plt.savefig("distance_first_2_days.png", dpi=200)

plt.show()

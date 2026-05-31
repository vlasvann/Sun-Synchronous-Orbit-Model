cmake_minimum_required(VERSION 3.16)

project(OrbitSimulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(orbit_core
    src/orbit/acceleration_model.cpp
    src/orbit/initial_conditions.cpp
    src/orbit/j2_perturbation.cpp
    src/orbit/orbit_simulation.cpp
    src/orbit/rk8_integrator.cpp
    src/orbit/orbital_elements.cpp
    src/orbit/sun_synchronous_orbit.cpp
)

target_include_directories(orbit_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_features(orbit_core PUBLIC cxx_std_17)

if(MSVC)
    target_compile_options(orbit_core PRIVATE /W4 /permissive-)
else()
    target_compile_options(orbit_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(orbit_simulation
    app/simulation_main.cpp
)

target_link_libraries(orbit_simulation
    PRIVATE
        orbit_core
)

# add_executable(orbit_visualization
#     app/visualization_main.cpp
# )

# target_link_libraries(orbit_visualization
#     PRIVATE
#         orbit_core
# )
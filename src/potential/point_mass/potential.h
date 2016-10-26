/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_POTENTIAL_POINT_MASS_H
#define SWIFT_POTENTIAL_POINT_MASS_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>

/* Local includes. */
#include "error.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

/**
 * @brief External Potential Properties - Point mass case
 */
struct external_potential {

  /*! Position of the point mass */
  double x, y, z;

  /*! Mass */
  double mass;

  /*! Time-step condition pre-factor */
  float timestep_mult;
};

/**
 * @brief Computes the time-step due to the acceleration from a point mass
 *
 * We pass in the time for simulations where the potential evolves with time.
 *
 * @param time The current time.
 * @param potential The properties of the externa potential.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static float external_gravity_timestep(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const,
    const struct gpart* restrict g) {

  const float G_newton = phys_const->const_newton_G;
  const float dx = g->x[0] - potential->x;
  const float dy = g->x[1] - potential->y;
  const float dz = g->x[2] - potential->z;
  const float rinv = 1.f / sqrtf(dx * dx + dy * dy + dz * dz);
  const float rinv2 = rinv * rinv;
  const float rinv3 = rinv2 * rinv;
  const float drdv = (g->x[0] - potential->x) * (g->v_full[0]) +
                     (g->x[1] - potential->y) * (g->v_full[1]) +
                     (g->x[2] - potential->z) * (g->v_full[2]);
  const float dota_x = G_newton * potential->mass * rinv3 *
                       (-g->v_full[0] + 3.f * rinv2 * drdv * dx);
  const float dota_y = G_newton * potential->mass * rinv3 *
                       (-g->v_full[1] + 3.f * rinv2 * drdv * dy);
  const float dota_z = G_newton * potential->mass * rinv3 *
                       (-g->v_full[2] + 3.f * rinv2 * drdv * dz);
  const float dota_2 = dota_x * dota_x + dota_y * dota_y + dota_z * dota_z;
  const float a_2 = g->a_grav[0] * g->a_grav[0] + g->a_grav[1] * g->a_grav[1] +
                    g->a_grav[2] * g->a_grav[2];

  return potential->timestep_mult * sqrtf(a_2 / dota_2);
}

/**
 * @brief Computes the gravitational acceleration of a particle due to a
 * point mass
 *
 * Note that the accelerations are multiplied by Newton's G constant later
 * on.
 *
 * We pass in the time for simulations where the potential evolves with time.
 *
 * @param time The current time.
 * @param potential The proerties of the external potential.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void external_gravity_acceleration(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const, struct gpart* restrict g) {

  const float dx = g->x[0] - potential->x;
  const float dy = g->x[1] - potential->y;
  const float dz = g->x[2] - potential->z;
  const float rinv = 1.f / sqrtf(dx * dx + dy * dy + dz * dz);
  const float rinv3 = rinv * rinv * rinv;

  g->a_grav[0] += -potential->mass * dx * rinv3;
  g->a_grav[1] += -potential->mass * dy * rinv3;
  g->a_grav[2] += -potential->mass * dz * rinv3;
}

/**
 * @brief Computes the gravitational potential energy of a particle in a point
 * mass potential.
 *
 * @param potential The #external_potential used in the run.
 * @param phys_const Physical constants in internal units.
 * @param g Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
external_gravity_get_potential_energy(
    const struct external_potential* potential,
    const struct phys_const* const phys_const, const struct gpart* g) {

  const float dx = g->x[0] - potential->x;
  const float dy = g->x[1] - potential->y;
  const float dz = g->x[2] - potential->z;
  const float rinv = 1. / sqrtf(dx * dx + dy * dy + dz * dz);
  return -phys_const->const_newton_G * potential->mass * r_inv;
}

/**
 * @brief Initialises the external potential properties in the internal system
 * of units.
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param s The #space we run in.
 * @param potential The external potential properties to initialize
 */
static INLINE void potential_init_backend(
    const struct swift_params* parameter_file,
    const struct phys_const* phys_const, const struct UnitSystem* us,
    const struct space* s, struct external_potential* potential) {

  potential->x =
      parser_get_param_double(parameter_file, "PointMassPotential:position_x");
  potential->y =
      parser_get_param_double(parameter_file, "PointMassPotential:position_y");
  potential->z =
      parser_get_param_double(parameter_file, "PointMassPotential:position_z");
  potential->mass =
      parser_get_param_double(parameter_file, "PointMassPotential:mass");
  potential->timestep_mult = parser_get_param_float(
      parameter_file, "PointMassPotential:timestep_mult");
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
static INLINE void potential_print_backend(
    const struct external_potential* potential) {

  message(
      "External potential is 'Point mass' with properties (x,y,z) = (%e, %e, "
      "%e), M = %e timestep multiplier = %e.",
      potential->x, potential->y, potential->z, potential->mass,
      potential->timestep_mult);
}

#endif /* SWIFT_POTENTIAL_POINT_MASS_H */

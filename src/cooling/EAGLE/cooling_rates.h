/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "interpolate.h"

/**
 * @brief calculates heating due to helium reionization in CGS units.
 *
 * @param z redshift
 * @param delta_z change in redshift over timestep
 * @param cooling the #cooling_function_data struct
 *
 * @return Cooling rate in CGS.
 */
__attribute__((always_inline)) INLINE double
eagle_helium_reionization_extraheat(double z, double delta_z,
                                    const struct cooling_function_data *restrict
                                        cooling) {
#ifdef SWIFT_DEBUG_CHECKS
  if (delta_z > 0.f) error("Invalid value for delta_z. Should be negative.");
#endif

  /* Recover the values we need */
  const double z_centre = cooling->He_reion_z_centre;
  const double z_sigma = cooling->He_reion_z_sigma;
  const double heat_cgs = cooling->He_reion_ev_pH;

  double extra_heat = 0.;

  /* Integral of the Gaussian between z and z - delta_z */
  extra_heat += erf((z - delta_z - z_centre) / (M_SQRT2 * z_sigma));
  extra_heat -= erf((z - z_centre) / (M_SQRT2 * z_sigma));

  /* Multiply by the normalisation factor */
  extra_heat *= heat_cgs * 0.5;

  return extra_heat;
}

/**
 * @brief Interpolates temperature from internal energy based on table and
 * calculates the size of the internal energy cell for the specified
 * internal energy. Returns log base 10 of temperature.
 *
 * @param log_10_u Log base 10 of internal energy.
 * @param redshift Current redshift.
 * @param compute_dT_du
 * @param dT_du Pointer to rate of change of Temperature with u.
 * @param z_i Redshift index
 * @param n_h_i Hydrogen number density index
 * @param He_i Helium fraction index
 * @param d_He Helium fraction offset
 * @param cooling #cooling_function_data structure
 */
__attribute__((always_inline)) INLINE double eagle_convert_u_to_temp(
    const double log_10_u_cgs, const float redshift, const int compute_dT_du,
    float *dT_du, int n_h_index, int He_index, float d_n_h, float d_He,
    const struct cooling_function_data *restrict cooling) {

  /* Get index of u along the table axis */
  int u_index;
  float d_u;
  get_index_1d(cooling->Therm, eagle_cooling_N_temperature, log_10_u_cgs,
               &u_index, &d_u);

  /* Interpolate temperature table to return temperature for current
   * internal energy (use 3D interpolation for high redshift table,
   * otherwise 4D) */
  float log_10_T;
  if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

    log_10_T = interpolation_3d(cooling->table.temperature,   /* */
                                n_h_index, He_index, u_index, /* */
                                d_n_h, d_He, d_u,             /* */
                                eagle_cooling_N_density,      /* */
                                eagle_cooling_N_He_frac,      /* */
                                eagle_cooling_N_temperature); /* */
  } else {

    log_10_T =
        interpolation_4d(cooling->table.temperature,                  /* */
                         /*z_index=*/0, n_h_index, He_index, u_index, /* */
                         cooling->dz, d_n_h, d_He, d_u,               /* */
                         eagle_cooling_N_loaded_redshifts,            /* */
                         eagle_cooling_N_density,                     /* */
                         eagle_cooling_N_He_frac,                     /* */
                         eagle_cooling_N_temperature);                /* */
  }

  if (compute_dT_du) {

    float log_10_T_high, log_10_T_low;

    /* Interpolate temperature table to return temperature for internal energy
     * at grid point above current internal energy for computing dT_du used for
     * calculation of dlambda_du in cooling.c (use 3D interpolation for high
     * redshift table, otherwise 4D) */
    if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

      log_10_T_high = interpolation_3d(cooling->table.temperature,   /* */
                                       n_h_index, He_index, u_index, /* */
                                       d_n_h, d_He, /*delta_u=*/1.f, /* */
                                       eagle_cooling_N_density,      /* */
                                       eagle_cooling_N_He_frac,      /* */
                                       eagle_cooling_N_temperature); /* */

    } else {

      log_10_T_high =
          interpolation_4d(cooling->table.temperature,                  /* */
                           /*z_index=*/0, n_h_index, He_index, u_index, /* */
                           cooling->dz, d_n_h, d_He, /*delta_u=*/1.f,   /* */
                           eagle_cooling_N_loaded_redshifts,            /* */
                           eagle_cooling_N_density,                     /* */
                           eagle_cooling_N_He_frac,                     /* */
                           eagle_cooling_N_temperature);                /* */
    }

    /* Interpolate temperature table to return temperature for internal energy
     * at grid point below current internal energy for computing dT_du used for
     * calculation of dlambda_du in cooling.c (use 3D interpolation for high
     * redshift table, otherwise 4D) */
    if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

      log_10_T_low = interpolation_3d(cooling->table.temperature,   /* */
                                      n_h_index, He_index, u_index, /* */
                                      d_n_h, d_He, /*delta_u=*/0.f, /* */
                                      eagle_cooling_N_density,      /* */
                                      eagle_cooling_N_He_frac,      /* */
                                      eagle_cooling_N_temperature); /* */

    } else {

      log_10_T_low =
          interpolation_4d(cooling->table.temperature,                  /* */
                           /*z_index=*/0, n_h_index, He_index, u_index, /* */
                           cooling->dz, d_n_h, d_He, /*delta_u=*/0.f,   /* */
                           eagle_cooling_N_loaded_redshifts,            /* */
                           eagle_cooling_N_density,                     /* */
                           eagle_cooling_N_He_frac,                     /* */
                           eagle_cooling_N_temperature);                /* */
    }

    /* Calculate dT/du */
    const float delta_u = exp(cooling->Therm[u_index + 1] * M_LN10) -
                          exp(cooling->Therm[u_index] * M_LN10);
    *dT_du =
        (exp(M_LN10 * log_10_T_high) - exp(M_LN10 * log_10_T_low)) / delta_u;
  }

  return log_10_T;
}

/**
 * @brief Calculates cooling rate for given internal energy by interpolating
 * EAGLE cooling tables.
 *
 * The tables depend on redshift, temperature, hydrogen number
 * density, helium fraction and metal abundance. Since only the
 * temperature changes when cooling a given particle, the redshift,
 * hydrogen number density and helium fraction indices and offsets
 * passed in.
 *
 * If the arguement dlambda_du is non-NULL, the routine also
 * calculates derivative of cooling rate with respect to internal
 * energy.
 *
 * If the argument element_lambda is non-NULL, the routine also
 * returns the cooling rate per element in the array.
 *
 * @param log10_u_cgs Log base 10 of internal energy per unit mass in CGS units.
 * @param redshift The current redshift
 * @param n_H_cgs The Hydrogen number density in CGS units.
 * @param solar_ratio Array of ratios of particle metal abundances
 * to solar metal abundances
 *
 * @param n_h_i Particle hydrogen number density index
 * @param d_n_h Particle hydrogen number density offset
 * @param He_i Particle helium fraction index
 * @param d_He Particle helium fraction offset
 * @param cooling Cooling data structure
 * @param phys_const Physical constants structure
 *
 * @param dlambda_du (return) Derivative of the cooling rate with respect to u.
 * @param element_lambda (return) Cooling rate from each element
 *
 * @return The cooling rate
 */
INLINE static double eagle_metal_cooling_rate(
    double log10_u_cgs, double redshift, double n_H_cgs,
    const float solar_ratio[chemistry_element_count + 2], int n_H_index,
    float d_n_h, int He_index, float d_He,
    const struct cooling_function_data *restrict cooling,
    const struct phys_const *phys_const, double *dlambda_du,
    double *element_lambda) {

  double h_plus_he_electron_abundance;

#ifdef TO_BE_DONE
  /* used for calculating dlambda_du */
  double temp_lambda_high = 0, temp_lambda_low = 0;
  double h_plus_he_electron_abundance_high = 0;
  double h_plus_he_electron_abundance_low = 0;
  double solar_electron_abundance_high = 0;
  double solar_electron_abundance_low = 0;
  double elem_cool_low = 0, elem_cool_high = 0;
#endif

  /* We only need dT_du if dLambda_du is non-NULL */
  const int compute_dT_du = (dlambda_du != NULL) ? 1 : 0;

  /* Temperature */
  float dT_du = -1.f;
  const double T =
      eagle_convert_u_to_temp(log10_u_cgs, redshift, compute_dT_du, &dT_du,
                              n_H_index, He_index, d_n_h, d_He, cooling);

  /* Get index along temperature dimension of the tables */
  int T_index;
  float d_T;
  get_index_1d(cooling->Temp, eagle_cooling_N_temperature, T, &T_index, &d_T);

#ifdef TO_BE_DONE
  /* Difference between entries on the temperature table around u */
  const float delta_T = exp(M_LN10 * cooling->Temp[T_index + 1]) -
                        exp(M_LN10 * cooling->Temp[T_index]);
#endif

  /* ------------------ */
  /* Metal-free cooling */
  /* ------------------ */

  double lambda_free = 0.;

  /* contribution to cooling and electron abundance from H, He. */
  if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

    /* If we're using the high redshift tables then we don't interpolate
     * in redshift */
    lambda_free = interpolation_3d(cooling->table.H_plus_He_heating, /* */
                                   n_H_index, He_index, T_index,     /* */
                                   d_n_h, d_He, d_T,                 /* */
                                   eagle_cooling_N_density,          /* */
                                   eagle_cooling_N_He_frac,          /* */
                                   eagle_cooling_N_temperature);     /* */

    h_plus_he_electron_abundance =
        interpolation_3d(cooling->table.H_plus_He_electron_abundance, /* */
                         n_H_index, He_index, T_index,                /* */
                         d_n_h, d_He, d_T,                            /* */
                         eagle_cooling_N_density,                     /* */
                         eagle_cooling_N_He_frac,                     /* */
                         eagle_cooling_N_temperature);                /* */

#ifdef TO_BE_DONE
    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du. Pass in NULL pointer for
     * dlambda_du in order to skip */
    if (dlambda_du != NULL) {
      temp_lambda_high = interpolation_3d(
          cooling->table.H_plus_He_heating, n_H_index, He_index, T_index, d_n_h,
          d_He, 1.f, cooling->N_nH, cooling->N_He, cooling->N_Temp);
      temp_lambda_low = interpolation_3d(
          cooling->table.H_plus_He_heating, n_H_index, He_index, T_index, d_n_h,
          d_He, 0.f, cooling->N_nH, cooling->N_He, cooling->N_Temp);
      h_plus_he_electron_abundance_high =
          interpolation_3d(cooling->table.H_plus_He_electron_abundance,
                           n_H_index, He_index, T_index, d_n_h, d_He, 1.f,
                           cooling->N_nH, cooling->N_He, cooling->N_Temp);
      h_plus_he_electron_abundance_low =
          interpolation_3d(cooling->table.H_plus_He_electron_abundance,
                           n_H_index, He_index, T_index, d_n_h, d_He, 0.f,
                           cooling->N_nH, cooling->N_He, cooling->N_Temp);
    }
#endif

  } else {

    /* Using normal tables, have to interpolate in redshift */
    lambda_free =
        interpolation_4d(cooling->table.H_plus_He_heating,            /* */
                         /*z_index=*/0, n_H_index, He_index, T_index, /* */
                         cooling->dz, d_n_h, d_He, d_T,               /* */
                         eagle_cooling_N_loaded_redshifts,            /* */
                         eagle_cooling_N_density,                     /* */
                         eagle_cooling_N_He_frac,                     /* */
                         eagle_cooling_N_temperature);                /* */

    h_plus_he_electron_abundance =
        interpolation_4d(cooling->table.H_plus_He_electron_abundance, /* */
                         /*z_index=*/0, n_H_index, He_index, T_index, /* */
                         cooling->dz, d_n_h, d_He, d_T,               /* */
                         eagle_cooling_N_loaded_redshifts,            /* */
                         eagle_cooling_N_density,                     /* */
                         eagle_cooling_N_He_frac,                     /* */
                         eagle_cooling_N_temperature);                /* */

#ifdef TO_BE_DONE
    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du */
    if (dlambda_du != NULL) {
      temp_lambda_high =
          interpolation_4d(cooling->table.H_plus_He_heating, 0, n_H_index,
                           He_index, T_index, cooling->dz, d_n_h, d_He, 1.f, 2,
                           cooling->N_nH, cooling->N_He, cooling->N_Temp);
      temp_lambda_low =
          interpolation_4d(cooling->table.H_plus_He_heating, 0, n_H_index,
                           He_index, T_index, cooling->dz, d_n_h, d_He, 0.f, 2,
                           cooling->N_nH, cooling->N_He, cooling->N_Temp);
      h_plus_he_electron_abundance_high = interpolation_4d(
          cooling->table.H_plus_He_electron_abundance, 0, n_H_index, He_index,
          T_index, cooling->dz, d_n_h, d_He, 1.f, 2, cooling->N_nH,
          cooling->N_He, cooling->N_Temp);
      h_plus_he_electron_abundance_low = interpolation_4d(
          cooling->table.H_plus_He_electron_abundance, 0, n_H_index, He_index,
          T_index, cooling->dz, d_n_h, d_He, 0.f, 2, cooling->N_nH,
          cooling->N_He, cooling->N_Temp);
    }
#endif
  }

#ifdef TO_BE_DONE
  if (dlambda_du != NULL) {
    *dlambda_du += (temp_lambda_high - temp_lambda_low) / delta_T * dT_du;
  }
#endif

  /* If we're testing cooling rate contributions write to array */
  if (element_lambda != NULL) {
    element_lambda[0] = lambda_free;
  }

  /* ------------------ */
  /* Compton cooling    */
  /* ------------------ */

  double lambda_Compton = 0.;

  /* inverse Compton cooling is not in collisional table
   * before reionisation so add now */

  if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1] ||
      redshift > cooling->reionisation_redshift) {

    const double zp1 = 1. + redshift;
    const double zp1p4 = zp1 * zp1 * zp1 * zp1;
    const double T_CMB = cooling->T_CMB_0 * zp1;

    lambda_Compton = -cooling->compton_rate_cgs * (T - T_CMB) * zp1p4 *
                     h_plus_he_electron_abundance / n_H_cgs;
  }

  if (element_lambda != NULL) {
    element_lambda[1] = lambda_Compton;
  }

  /* --------------------------------- */
  /* Electron abundance ratio to solar */
  /* --------------------------------- */

  double solar_electron_abundance;

  if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

    /* If we're using the high redshift tables then we don't interpolate
     * in redshift */
    solar_electron_abundance =
        interpolation_2d(cooling->table.electron_abundance, /* */
                         n_H_index, T_index,                /* */
                         d_n_h, d_T,                        /* */
                         eagle_cooling_N_density,           /* */
                         eagle_cooling_N_temperature);      /* */

#ifdef TO_BE_DONE
    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du */
    if (dlambda_du != NULL) {
      solar_electron_abundance_high =
          interpolation_2d(cooling->table.electron_abundance, n_H_index,
                           T_index, d_n_h, 1.f, cooling->N_nH, cooling->N_Temp);
      solar_electron_abundance_low =
          interpolation_2d(cooling->table.electron_abundance, n_H_index,
                           T_index, d_n_h, 0.f, cooling->N_nH, cooling->N_Temp);
    }
#endif

  } else {

    /* Using normal tables, have to interpolate in redshift */
    solar_electron_abundance =
        interpolation_3d(cooling->table.electron_abundance, /* */
                         /*z_index=*/0, n_H_index, T_index, /* */
                         cooling->dz, d_n_h, d_T,           /* */
                         eagle_cooling_N_loaded_redshifts,  /* */
                         eagle_cooling_N_density,           /* */
                         eagle_cooling_N_temperature);      /* */

#ifdef TO_BE_DONE
    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du */
    if (dlambda_du != NULL) {
      solar_electron_abundance_high = interpolation_3d(
          cooling->table.electron_abundance, 0, n_H_index, T_index, cooling->dz,
          d_n_h, 1.f, 2, cooling->N_nH, cooling->N_Temp);
      solar_electron_abundance_low = interpolation_3d(
          cooling->table.electron_abundance, 0, n_H_index, T_index, cooling->dz,
          d_n_h, 0.f, 2, cooling->N_nH, cooling->N_Temp);
    }
#endif
  }

  /* ------------- */
  /* Metal cooling */
  /* ------------- */

  /* for each element the cooling rate is multiplied by the ratio of H, He
   * electron abundance to solar electron abundance then by the ratio of the
   * particle metal abundance to solar metal abundance. */

  const double abundance_ratio =
      h_plus_he_electron_abundance / solar_electron_abundance;

  double lambda_metal[eagle_cooling_N_metal];

  if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

    for (int elem = 0; elem < eagle_cooling_N_metal; elem++) {

      lambda_metal[elem] =
          interpolation_3d_no_x(cooling->table.metal_heating,   /* */
                                elem, n_H_index, T_index,       /* */
                                /*delta_elem=*/0.f, d_n_h, d_T, /* */
                                eagle_cooling_N_metal,          /* */
                                eagle_cooling_N_density,        /* */
                                eagle_cooling_N_temperature);   /* */

      lambda_metal[elem] *= abundance_ratio;
      lambda_metal[elem] *= solar_ratio[elem + 2];

#ifdef TO_BE_DONE
      /* compute values at temperature gridpoints above and below input
       * temperature for calculation of dlambda_du */
      if (dlambda_du != NULL) {
        elem_cool_high = interpolation_3d_no_x(
            cooling->table.metal_heating, elem, n_H_index, T_index, 0.f, d_n_h,
            1.f, cooling->N_Elements, cooling->N_nH, cooling->N_Temp);

        elem_cool_low = interpolation_3d_no_x(
            cooling->table.metal_heating, elem, n_H_index, T_index, 0.f, d_n_h,
            0.f, cooling->N_nH, cooling->N_Temp, cooling->N_Elements);

        *dlambda_du += (elem_cool_high * h_plus_he_electron_abundance_high /
                            solar_electron_abundance_high -
                        elem_cool_low * h_plus_he_electron_abundance_low /
                            solar_electron_abundance_low) /
                       delta_T * dT_du * solar_ratio[elem + 2];
      }
#endif
    }

  } else {

    for (int elem = 0; elem < eagle_cooling_N_metal; elem++) {

      lambda_metal[elem] = interpolation_4d_no_x(
          cooling->table.metal_heating,                /* */
          elem, /*z_index=*/0, n_H_index, T_index,     /* */
          /*delta_elem=*/0.f, cooling->dz, d_n_h, d_T, /* */
          eagle_cooling_N_metal,                       /* */
          eagle_cooling_N_loaded_redshifts,            /* */
          eagle_cooling_N_density,                     /* */
          eagle_cooling_N_temperature);                /* */

      lambda_metal[elem] *= abundance_ratio;
      lambda_metal[elem] *= solar_ratio[elem + 2];

#ifdef TO_BE_DONE
      /* compute values at temperature gridpoints above and below input
       * temperature for calculation of dlambda_du */
      if (dlambda_du != NULL) {
        elem_cool_high = interpolation_4d_no_x(
            cooling->table.metal_heating, elem, 0, n_H_index, T_index, 0.,
            cooling->dz, d_n_h, 1.f, cooling->N_Elements, 2, cooling->N_nH,
            cooling->N_Temp);

        elem_cool_low = interpolation_4d_no_x(
            cooling->table.metal_heating, elem, 0, n_H_index, T_index, 0.,
            cooling->dz, d_n_h, 0.f, cooling->N_Elements, 2, cooling->N_nH,
            cooling->N_Temp);

        *dlambda_du += (elem_cool_high * h_plus_he_electron_abundance_high /
                            solar_electron_abundance_high -
                        elem_cool_low * h_plus_he_electron_abundance_low /
                            solar_electron_abundance_low) /
                       delta_T * dT_du * solar_ratio[elem + 2];
      }
#endif
    }
  }

  if (element_lambda != NULL) {
    for (int elem = 0; elem < eagle_cooling_N_metal; ++elem) {
      element_lambda[elem + 2] = lambda_metal[elem];
    }
  }

  /* Sum up all the contributions */
  double cooling_rate = lambda_free + lambda_Compton;
  for (int elem = 0; elem < eagle_cooling_N_metal; ++elem) {
    cooling_rate += lambda_metal[elem];
  }

  return cooling_rate;
}

/**
 * @brief Wrapper function used to calculate cooling rate and dLambda_du.
 * Table indices and offsets for redshift, hydrogen number density and
 * helium fraction are passed it so as to compute them only once per particle.
 *
 * @param log_u_cgs Natural log of internal energy per unit mass in CGS units.
 * @param redshift The current redshift.
 * @param n_H_cgs Hydrogen number density in CGS units.
 * @param abundance_ratio Ratio of element abundance to solar.
 *
 * @param n_h_i Particle hydrogen number density index
 * @param d_n_h Particle hydrogen number density offset
 * @param He_i Particle helium fraction index
 * @param d_He Particle helium fraction offset
 * @param cooling #cooling_function_data structure
 * @param phys_const #phys_const structure
 *
 * @param dLambdaNet_du (return) Derivative of the cooling rate with respect to
 * u.
 *
 * @return The cooling rate
 */
INLINE static double eagle_cooling_rate(
    double log_u_cgs, double redshift, double n_H_cgs,
    const float abundance_ratio[chemistry_element_count + 2], int n_H_index,
    float d_n_h, int He_index, float d_He,
    const struct cooling_function_data *restrict cooling,
    const struct phys_const *phys_const, double *dLambdaNet_du) {

  return eagle_metal_cooling_rate(log_u_cgs / M_LN10, redshift, n_H_cgs,
                                  abundance_ratio, n_H_index, d_n_h, He_index,
                                  d_He, cooling, phys_const, dLambdaNet_du,
                                  /*element_lambda=*/NULL);
}

/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *                    Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_GEAR_STARS_IACT_H
#define SWIFT_GEAR_STARS_IACT_H

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sparticle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param xp Extra particle data
 * @param ti_current Current integer time value
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_stars_density(
    float r2, const float *dx, float hi, float hj, struct spart *restrict si,
    const struct part *restrict pj, const struct cosmology *restrict cosmo,
    const struct stars_props *restrict stars_properties,
    struct xpart *restrict xp, integertime_t ti_current) {

  float wi, wi_dx;

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the number of neighbours */
  si->density.wcount += wi;
  si->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  /* Compute contribution to the density */
  si->density.rho_gas += pj->mass * wi;


#ifdef DEBUG_INTERACTIONS_STARS
  /* Update ngb counters */
  if (si->num_ngb_density < MAX_NUM_OF_NEIGHBOURS_STARS)
    si->ids_ngbs_density[si->num_ngb_density] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_density;
#endif
}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 * Used for updating properties of gas particles neighbouring a star particle
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle.
 * @param pj Second (gas) particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param xp Extra particle data
 * @param ti_current Current integer time used value for seeding random number
 * generator
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_stars_feedback(
    float r2, const float *dx, float hi, float hj, struct spart *restrict si,
    struct part *restrict pj, const struct cosmology *restrict cosmo,
    const struct stars_props *restrict stars_properties,
    struct xpart *restrict xp, integertime_t ti_current) {

  const float mj = hydro_get_mass(pj);
  const float rhoj = hydro_get_comoving_density(pj);
  const float r = sqrtf(r2);
  const float ri = 1.f / r;

  /* Get the kernel for hi. */
  float hi_inv = 1.0f / hi;
  float hi_inv_dim = pow_dimension(hi_inv);       /* 1/h^d */
  float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  float wi_dr = hid_inv * wi_dx;

  /* Compute dv dot r */
  float dvdr = (si->v[0] - pj->v[0]) * dx[0] + (si->v[1] - pj->v[1]) * dx[1] +
               (si->v[2] - pj->v[2]) * dx[2];

  /* Get the time derivative for h. */
  si->feedback.h_dt -= mj * dvdr * ri / rhoj * wi_dr;

#ifdef DEBUG_INTERACTIONS_STARS
  /* Update ngb counters */
  if (si->num_ngb_force < MAX_NUM_OF_NEIGHBOURS_STARS)
    si->ids_ngbs_force[si->num_ngb_force] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_force;
#endif

  /* Does the star need to explode? */
  if (si->has_exploded > 0)
    return;

  si->has_exploded = -1;

  /* Get a few variables */
  const double e_sn = stars_properties->feedback.energy_per_supernovae;
  const double m_ej = stars_properties->feedback.mass_ejected;
  const double e_sn_normalized = e_sn * m_ej / si->density.rho_gas;
  const double u_old = hydro_get_physical_internal_energy(pj, xp, cosmo);

  /* Mass received */
  const double wij = mj * wi * hi_inv_dim / si->density.rho_gas;
  const double dm = m_ej * wij;

  /* Energy received */
  const double du = m_ej * e_sn_normalized * wij / (mj + dm);

  /* TODO updates everything in between steps */
  hydro_set_mass(pj, mj + dm);
  stars_set_mass(si, si->mass - dm);
  /* error("%g %g %g %g", u_old, du, si->density.rho_gas, e_sn); */
  /* hydro_set_physical_internal_energy(pj, xp, cosmo, u_old + du); */
  /* hydro_set_drifted_physical_internal_energy(pj, cosmo, u_old + du); */

  /* // TODO activate particle */

  /* /\* Compute the norm of the speed. *\/ */
  /* float vi2 = 0.; */
  /* float vj2 = 0.; df*/
  /* for(int i = 0; i < 3; i++) { */
  /*   vi2 += si->v[i] * si->v[i]; */
  /*   vj2 += pj->v[i] * pj->v[i]; */
  /* } */
  
  /* float fac = mj / (mj + dm) + wi * h_inv_dim * */
  /*   m_ej * vi2 * vj2 / (mj + dm); */
  /* fac = 1 - sqrt(fac); */
  /* for(int i = 0; i < 3; i++) { */
  /*   // TODO */
  /* } */

}


#endif // SWIFT_GEAR_STARS_IACT_H

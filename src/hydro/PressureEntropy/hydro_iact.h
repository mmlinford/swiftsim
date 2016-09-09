/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_PRESSURE_ENTROPY_HYDRO_IACT_H
#define SWIFT_PRESSURE_ENTROPY_HYDRO_IACT_H

/**
 * @file PressureEntropy/hydro_iact.h
 * @brief Pressure-Entropy implementation of SPH (Neighbour loop equations)
 *
 * The thermal variable is the entropy (S) and the entropy is smoothed over
 * contact discontinuities to prevent spurious surface tension.
 *
 * Follows Hopkins, P., MNRAS, 2013, Volume 428, Issue 4, pp. 2840-2856
 */

/**
 * @brief Density loop
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wi_dx;
  float wj, wj_dx;
  /* float dv[3], curlvr[3]; */

  /* Get the masses. */
  const float mi = pi->mass;
  const float mj = pj->mass;

  /* Get the entropies. */
  const float entropy_i = pi->entropy;
  const float entropy_j = pj->entropy;

  /* Get r and r inverse. */
  const float r = sqrtf(r2);
  /* const float r_inv = 1.0f / r; */

  /* Compute the kernel function for pi */
  const float hi_inv = 1.f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the density */
  pi->rho += mj * wi;

  /* Weighted pressure */
  pi->weightedPressure += mj * pow_one_over_gamma(pj->entropy) * wi;
  pi->density.weightedPressure_dh -=
      mj * pow_one_over_gamma(entropy_j) * (hydro_dimension * wi + ui * wi_dx);

  /* Compute contribution to the number of neighbours */
  pi->density.wcount += wi;
  pi->density.wcount_dh -= ui * wi_dx;  //(hydro_dimension * wi + ui * wi_dx);

  /* Compute the kernel function for pj */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  /* Compute contribution to the density */
  pj->rho += mi * wj;
  // pj->rho_dh -= mi * (hydro_dimension * wj + uj * wj_dx);

  /* Weighted pressure */
  pj->weightedPressure += mi * pow_one_over_gamma(pi->entropy) * wj;
  pj->density.weightedPressure_dh -=
      mi * pow_one_over_gamma(entropy_i) * (hydro_dimension * wj + uj * wj_dx);

  /* Compute contribution to the number of neighbours */
  pj->density.wcount += wj;
  pj->density.wcount_dh -= uj * wj_dx;  //(hydro_dimension * wj + uj * wj_dx);

  /* const float faci = mj * wi_dx * r_inv; */
  /* const float facj = mi * wj_dx * r_inv; */

  /* /\* Compute dv dot r *\/ */
  /* dv[0] = pi->v[0] - pj->v[0]; */
  /* dv[1] = pi->v[1] - pj->v[1]; */
  /* dv[2] = pi->v[2] - pj->v[2]; */
  /* const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2]; */

  /* pi->density.div_v -= faci * dvdr; */
  /* pj->density.div_v -= facj * dvdr; */

  /* /\* Compute dv cross r *\/ */
  /* curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1]; */
  /* curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2]; */
  /* curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0]; */

  /* pi->density.rot_v[0] += faci * curlvr[0]; */
  /* pi->density.rot_v[1] += faci * curlvr[1]; */
  /* pi->density.rot_v[2] += faci * curlvr[2]; */

  /* pj->density.rot_v[0] += facj * curlvr[0]; */
  /* pj->density.rot_v[1] += facj * curlvr[1]; */
  /* pj->density.rot_v[2] += facj * curlvr[2]; */
}

/**
 * @brief Density loop (Vectorized version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_vec_density(
    float *R2, float *Dx, float *Hi, float *Hj, struct part **pi,
    struct part **pj) {

  error("Vectorized version of Pressure-Entropy SPH routine not existant yet.");
}

/**
 * @brief Density loop (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wi_dx;
  /* float dv[3], curlvr[3]; */

  /* Get the masses. */
  const float mj = pj->mass;

  /* Get the entropies. */
  const float entropy_j = pj->entropy;

  /* Get r and r inverse. */
  const float r = sqrtf(r2);
  /* const float ri = 1.0f / r; */

  /* Compute the kernel function */
  const float h_inv = 1.0f / hi;
  const float u = r * h_inv;
  kernel_deval(u, &wi, &wi_dx);

  /* Compute contribution to the density */
  pi->rho += mj * wi;

  /* Weighted pressure */
  pi->weightedPressure += mj * pow_one_over_gamma(pj->entropy) * wi;
  pi->density.weightedPressure_dh -=
      mj * pow_one_over_gamma(entropy_j) * (hydro_dimension * wi + u * wi_dx);

  /* Compute contribution to the number of neighbours */
  pi->density.wcount += wi;
  pi->density.wcount_dh -= u * wi_dx;  //(hydro_dimension * wi + u * wi_dx);

  /* const float fac = mj * wi_dx * ri; */

  /* /\* Compute dv dot r *\/ */
  /* dv[0] = pi->v[0] - pj->v[0]; */
  /* dv[1] = pi->v[1] - pj->v[1]; */
  /* dv[2] = pi->v[2] - pj->v[2]; */
  /* const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2]; */
  /* pi->density.div_v -= fac * dvdr; */

  /* /\* Compute dv cross r *\/ */
  /* curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1]; */
  /* curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2]; */
  /* curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0]; */

  /* pi->density.rot_v[0] += fac * curlvr[0]; */
  /* pi->density.rot_v[1] += fac * curlvr[1]; */
  /* pi->density.rot_v[2] += fac * curlvr[2]; */
}

/**
 * @brief Density loop (non-symmetric vectorized version)
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_vec_density(float *R2, float *Dx, float *Hi, float *Hj,
                               struct part **pi, struct part **pj) {

  error("Vectorized version of Pressure-Entropy SPH routine not existant yet.");
}

/**
 * @brief Force loop
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wj, wi_dx, wj_dx;

  const float fac_mu = 1.f; /* Will change with cosmological integration */

  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  /* Get some values in local variables. */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float entropy_i = pi->entropy;
  const float entropy_j = pj->entropy;
  const float entropy_product = pow_one_over_gamma(entropy_i * entropy_j);

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Compute gradient terms */
  const float f_ij =
      1.f;  // - pi->force.f_ij / (mj * pow_one_over_gamma(entropy_j));
  const float f_ji =
      1.f;  // - pj->force.f_ij / (mi * pow_one_over_gamma(entropy_i));

  /* Compute sound speeds */
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Balsara term */
  /* const float balsara_i = pi->force.balsara; */
  /* const float balsara_j = pj->force.balsara; */

  /* Are the particles moving towards each others ? */
  const float omega_ij = fminf(dvdr, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float v_sig = ci + cj - 3.f * mu_ij;

  /* Now construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.5f * const_viscosity_alpha * v_sig * mu_ij / rho_ij;

  /* Now, convolve with the kernel */
  const float visc_term = 0.5f * visc * (wi_dr + wj_dr) * r_inv;
  const float sph_term = entropy_product *
                         (f_ij * pi->force.pressure_term * wi_dr +
                          f_ji * pj->force.pressure_term * wj_dr) *
                         r_inv;

  /* if(pi->entropy != pj->entropy) { */
  /* message("Si = %f Pi = %f, Pi_term = %f", pi->entropy, pi->weightedPressure,
   * pi->force.pressure_term); */
  /* message("Sj = %f Pj = %f, Pj_term = %f", pj->entropy, pj->weightedPressure,
   * pj->force.pressure_term); */
  /* error("done"); */
  /* } */

  /* Eventually got the acceleration */
  const float acc = visc_term + sph_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  pj->a_hydro[0] += mi * acc * dx[0];
  pj->a_hydro[1] += mi * acc * dx[1];
  pj->a_hydro[2] += mi * acc * dx[2];

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;
  pj->force.h_dt -= mi * dvdr * r_inv / rhoi * wj_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = fmaxf(pi->force.v_sig, v_sig);
  pj->force.v_sig = fmaxf(pj->force.v_sig, v_sig);

  /* Change in entropy */
  pi->entropy_dt += mj * visc_term * dvdr;
  pj->entropy_dt += mi * visc_term * dvdr;
}

/**
 * @brief Force loop (Vectorized version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_vec_force(
    float *R2, float *Dx, float *Hi, float *Hj, struct part **pi,
    struct part **pj) {

  error("Vectorized version of Pressure-Entropy SPH routine not existant yet.");
}

/**
 * @brief Force loop (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wj, wi_dx, wj_dx;

  const float fac_mu = 1.f; /* Will change with cosmological integration */

  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  /* Get some values in local variables. */
  // const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float entropy_i = pi->entropy;
  const float entropy_j = pj->entropy;
  const float entropy_product = pow_one_over_gamma(entropy_i * entropy_j);

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Compute gradient terms */
  const float f_ij =
      1.f;  // - pi->force.f_ij / (mj * pow_one_over_gamma(entropy_j));
  const float f_ji =
      1.f;  // - pj->force.f_ij / (mi * pow_one_over_gamma(entropy_i));

  /* Compute sound speeds */
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Balsara term */
  /* const float balsara_i = pi->force.balsara; */
  /* const float balsara_j = pj->force.balsara; */

  /* Are the particles moving towards each others ? */
  const float omega_ij = fminf(dvdr, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float v_sig = ci + cj - 3.f * mu_ij;

  /* Now construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.5f * const_viscosity_alpha * v_sig * mu_ij / rho_ij;

  /* Now, convolve with the kernel */
  const float visc_term = 0.5f * visc * (wi_dr + wj_dr) * r_inv;
  const float sph_term = entropy_product *
                         (f_ij * pi->force.pressure_term * wi_dr +
                          f_ji * pj->force.pressure_term * wj_dr) *
                         r_inv;

  /* Eventually got the acceleration */
  const float acc = visc_term + sph_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = fmaxf(pi->force.v_sig, v_sig);

  /* Change in entropy */
  pi->entropy_dt += mj * visc_term * dvdr;
}

/**
 * @brief Force loop (Vectorized non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_vec_force(
    float *R2, float *Dx, float *Hi, float *Hj, struct part **pi,
    struct part **pj) {

  error("Vectorized version of Pressure-Entropy SPH routine not existant yet.");
}

#endif /* SWIFT_PRESSURE_ENTROPY_HYDRO_IACT_H */

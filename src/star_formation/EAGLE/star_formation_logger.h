/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
 *******************************************************************************/
#ifndef SWIFT_SCHAYE_STARFORMATION_LOGGER_H
#define SWIFT_SCHAYE_STARFORMATION_LOGGER_H

/* Some standard headers */
#include <stdlib.h>

/* Local includes */
#include "cell.h"
#include "hydro.h"
#include "part.h"
#include "star_formation_logger_struct.h"
#include "units.h"

/**
 * @brief Update the stellar mass in the current cell after creating
 * the new star particle spart sp
 *
 * @param sp new created star particle
 * @param sf the star_formation_history struct of the current cell
 */
INLINE static void star_formation_logger_log_new_spart(
    const struct spart *sp, struct star_formation_history *sf) {

  /* Add mass of created sparticle to the total stellar mass in this cell*/
  sf->new_stellar_mass += sp->mass;
}

/**
 * @brief Initialize the star formation history struct if the cell is active
 *
 * @param sf the star_formation_history struct we want to initialize
 */
INLINE static void star_formation_logger_log_active_cell(
    struct star_formation_history *sf) {

  /* Initialize the stellar mass to zero*/
  sf->new_stellar_mass = 0.f;

  /* Initialize the active SFR */
  sf->SFR_active = 0.f;

  /* Initialize the SFR*dt active */
  sf->SFRdt_active = 0.f;

  /* Initialize the inactive SFR */
  sf->SFR_inactive = 0.f;
}

/**
 * @brief Initialize the star formation history struct in the case the cell is
 * inactive
 *
 * @param sf the star_formation_history struct we want to initialize
 */
INLINE static void star_formation_logger_log_inactive_cell(
    struct star_formation_history *sf) {

  /* Initialize the stellar mass to zero*/
  sf->new_stellar_mass = 0.f;

  /* The active SFR becomes the inactive SFR */
  sf->SFR_inactive += sf->SFR_active;

  /* Initialize the active SFR */
  sf->SFR_active = 0.f;

  /* Initialize the SFR*dt active */
  sf->SFRdt_active = 0.f;
}

/**
 * @brief function to add the progeny SFH to the parent SFH.
 *
 * @param sf parent SFH struct
 * @param sfprogeny progeny SFH struct
 */
INLINE static void star_formation_logger_log_progeny_cell(
    struct star_formation_history *sf,
    const struct star_formation_history *sfprogeny) {

  /* Add the new stellar mass from the progeny */
  sf->new_stellar_mass += sfprogeny->new_stellar_mass;

  /* Add the SFR of the progeny */
  sf->SFR_active += sfprogeny->SFR_active;

  /* Add the SFR * dt of the progeny */
  sf->SFRdt_active += sfprogeny->SFRdt_active;

  /* Add the Inactive SFR of the progeny */
  sf->SFR_inactive += sfprogeny->SFR_inactive;
}

/**
 * @brief Get the total star formation in this cell and add it to the star
 * formation history struct in the #engine
 *
 * @param c the cell of which we want to know the star formation
 * @param sf the star formation structure to which we want to add the star
 * formation
 */
INLINE static void star_formation_logger_add(
    struct cell *c, struct star_formation_history *sf) {

  /* Get the star formation history struct from the cell */
  struct star_formation_history *sfcell = &c->stars.sfh;

  /* Update the SFH structure */
  sf->new_stellar_mass += sfcell->new_stellar_mass;

  sf->SFR_active += sfcell->SFR_active;

  sf->SFRdt_active += sfcell->SFRdt_active;

  sf->SFR_inactive += sfcell->SFR_inactive;
}

/**
 * @brief add the star formation to the parent cell in the #engine
 *
 * @param c the cell for which we want to add the star formation
 * @param sf the combined star formation history of the progeny
 */
INLINE static void star_formation_logger_assign(
    struct cell *c, const struct star_formation_history *sf) {

  /* Get the star formation history from the cell */
  struct star_formation_history *sfcell = &c->stars.sfh;

  /* Update the SFH structure */
  sfcell->new_stellar_mass += sf->new_stellar_mass;

  sfcell->SFR_active += sf->SFR_active;

  sfcell->SFRdt_active += sf->SFR_active;

  sfcell->SFR_inactive += sf->SFR_inactive;
}

/**
 * @brief add the star formation to the parent cell in the #engine
 *
 * @param c the cell for which we want to add the star formation
 * @param sf the combined star formation history of the progeny
 */
INLINE static void star_formation_logger_assign2(
    struct cell *c, const struct star_formation_history *sf) {

  /* Get the star formation history from the cell */
  struct star_formation_history *sfcell = &c->stars.sfh;

  /* Update the SFH structure */
  sfcell->new_stellar_mass = sf->new_stellar_mass;

  sfcell->SFR_active = sf->SFR_active;

  sfcell->SFRdt_active = sf->SFR_active;

  sfcell->SFR_inactive = sf->SFR_inactive;
}

/**
 * @brief Initialize the star formation history structure in the #engine
 *
 * @param The pointer to the star formation history structure
 */
INLINE static void star_formation_logger_init_engine(
    struct star_formation_history *sfh) {

  /* Initialize the collecting SFH structure to zero */
  sfh->new_stellar_mass = 0.f;

  sfh->SFR_active = 0.f;

  sfh->SFRdt_active = 0.f;

  sfh->SFR_inactive = 0.f;
}

/**
 * @brief Write the final SFH to a file
 *
 * @param fp the file pointer
 * @param time the simulation time
 * @param a the scale factor
 * @param z the redshift
 * @param sf the star_formation_history struct
 */
INLINE static void star_formation_logger_write_to_log_file(
    FILE *fp, const double time, const double a, const double z,
    const struct star_formation_history sf, const int step, const double total_SFH, const double active_SFH, const double inactive_SFH) {

  /* Calculate the total SFR */
  const float totalSFR = sf.SFR_active + sf.SFR_inactive;
  fprintf(fp, "%6d %16e %12.7f %12.7f %14e %14e %14e %14e %14e %14e %14e\n",step, time, a, z,
          sf.new_stellar_mass, sf.SFR_active, sf.SFRdt_active, totalSFR, total_SFH, active_SFH, inactive_SFH);
}

/**
 * @brief Initialize the SFH logger file
 *
 * @param fp the file pointer
 * @param us The current internal system of units.
 * @param phys_const Physical constants in internal units
 */
INLINE static void star_formation_logger_init_log_file(
    FILE *fp, const struct unit_system *restrict us,
    const struct phys_const *phys_const) {

  /* Write some general text to the logger file */
  fprintf(fp, "# Star Formation History Logger file\n");
  fprintf(fp, "######################################################\n");
  fprintf(fp, "# The quantities are all given in internal physical units!\n");
  fprintf(fp, "#\n");
  fprintf(fp, "# (0) Simulation step\n");
  fprintf(fp, "# (1) Time since Big Bang (cosmological run), Time since start of the simulation (non-cosmological run).\n");
  fprintf(fp, "#     Unit = %e seconds\n", us->UnitTime_in_cgs);
  fprintf(fp, "#     Unit = %e yr or %e Myr\n", 1.f / phys_const->const_year,
          1.f / phys_const->const_year / 1e6);
  fprintf(fp, "# (2) Scale factor (no unit)\n");
  fprintf(fp, "# (3) Redshift     (no unit)\n");
  fprintf(fp, "# (4) Total mass stars formed in the current time-step.\n");
  fprintf(fp, "#     Unit = %e gram\n", us->UnitMass_in_cgs);
  fprintf(fp, "#     Unit = %e solar mass\n",
          1.f / phys_const->const_solar_mass);
  fprintf(fp, "# (5) The total SFR of all the active particles.\n");
  fprintf(fp, "#     Unit = %e gram/s\n",
          us->UnitMass_in_cgs / us->UnitTime_in_cgs);
  fprintf(fp, "#     Unit = %e Msol/yr\n",
          phys_const->const_year / phys_const->const_solar_mass);
  fprintf(
      fp,
      "# (6) The star formation rate (SFR) of active particles multiplied by their time-step size.\n");
  fprintf(fp, "#     Unit = %e gram\n", us->UnitMass_in_cgs);
  fprintf(fp, "#     Unit = %e solar mass\n",
          1.f / phys_const->const_solar_mass);
  fprintf(fp, "# (7) The total SFR of all the particles in the simulation.\n");
  fprintf(fp, "#     Unit = %e gram/s\n",
          us->UnitMass_in_cgs / us->UnitTime_in_cgs);
  fprintf(fp, "#     Unit = %e Msol/yr\n",
          phys_const->const_year / phys_const->const_solar_mass);
  fprintf(fp, "#\n");
  fprintf(fp,
          "# (0)         (1)            (2)          (3)            (4)           "
          " (5)            (6)            (7)\n");
  fprintf(fp,
          "#            Time             a            z        total M_stars  SFR "
          "(active)  SFR*dt (active)  SFR (total)\n");
}

/**
 * @brief Add the SFR tracer to the total active SFR of this cell
 *
 * @param p the #part
 * @param xp the #xpart
 * @param sf the SFH logger struct
 */
INLINE static void star_formation_logger_log_active_part(
    const struct part *p, const struct xpart *xp,
    struct star_formation_history *sf, const double dt_star) {

  /* Add the SFR to the logger file */
  sf->SFR_active += xp->sf_data.SFR;

  /* Update the active SFR*dt */
  sf->SFRdt_active += xp->sf_data.SFR * dt_star;
}

/**
 * @brief Add the SFR tracer to the total inactive SFR of this cell as long as
 * the SFR tracer is larger than 0
 *
 * @param p the #part
 * @param xp the #xpart
 * @param sf the SFH logger struct
 */
INLINE static void star_formation_logger_log_inactive_part(
    const struct part *p, const struct xpart *xp,
    struct star_formation_history *sf) {

  /* Add the SFR to the logger file */
  sf->SFR_inactive += max(xp->sf_data.SFR, 0.f);
}

INLINE static void star_formation_logger_clean_split_cell(
    struct star_formation_history *sf){
  sf->new_stellar_mass = 0.f;
  sf->SFR_active = 0.f;
  sf->SFRdt_active = 0.f;
  sf->SFR_inactive = 0.f;
}

#endif /* SWIFT_SCHAYE_STARFORMATION_LOGGER_H */

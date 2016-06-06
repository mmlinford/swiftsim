/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

__attribute__((always_inline)) INLINE static void hydro_debug_particle(
    struct part* p, struct xpart* xp) {
  printf(
      "x=[%.3e,%.3e,%.3e], "
      "v=[%.3e,%.3e,%.3e],v_full=[%.3e,%.3e,%.3e] \n a=[%.3e,%.3e,%.3e], "
      "u_full=%.3e, u=%.3e, du/dt=%.3e v_sig=%.3e, P=%.3e\n"
      "h=%.3e, dh/dt=%.3e "
      "wcount=%d, m=%.3e, dh_drho=%.3e, rho=%.3e, t_begin=%d, t_end=%d\n",
      p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2], xp->v_full[0],
      xp->v_full[1], xp->v_full[2], p->a_hydro[0], p->a_hydro[1], p->a_hydro[2],
      xp->u_full, p->u, p->u_dt, p->force.v_sig, p->force.pressure, p->h,
      p->h_dt, (int)p->density.wcount, p->mass, p->rho_dh, p->rho, p->ti_begin,
      p->ti_end);
}

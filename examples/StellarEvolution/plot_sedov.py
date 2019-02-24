###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 #                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 # 
 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU Lesser General Public License as published
 # by the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 # 
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 # 
 # You should have received a copy of the GNU Lesser General Public License
 # along with this program.  If not, see <http://www.gnu.org/licenses/>.
 # 
 ##############################################################################

# Computes the analytical solution of the 3D Sedov blast wave.
# The script works for a given initial box and dumped energy and computes the solution at a later time t.

# Parameters
rho_0 = 1.          # Background Density
P_0 = 1.e-6         # Background Pressure
E_0 = 1.            # Energy of the explosion
gas_gamma = 5./3.   # Gas polytropic index


# ---------------------------------------------------------------
# Don't touch anything after this.
# ---------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
from pylab import *
from scipy import stats
import h5py

# Plot parameters
params = {'axes.labelsize': 10,
'axes.titlesize': 10,
'font.size': 12,
'legend.fontsize': 12,
'xtick.labelsize': 10,
'ytick.labelsize': 10,
'text.usetex': True,
 'figure.figsize' : (9.90,6.45),
'figure.subplot.left'    : 0.1,
'figure.subplot.right'   : 0.99,
'figure.subplot.bottom'  : 0.1,
'figure.subplot.top'     : 0.95,
'figure.subplot.wspace'  : 0.2,
'figure.subplot.hspace'  : 0.2,
'lines.markersize' : 6,
'lines.linewidth' : 3.,
'text.latex.unicode': True
}
rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Times']})


#snap = int(sys.argv[1])
snapshot = sys.argv[1]


# Read the simulation data
sim = h5py.File(snapshot, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]

# Units
unit_length_in_cgs = sim["/Units"].attrs["Unit length in cgs (U_L)"]
unit_mass_in_cgs = sim["/Units"].attrs["Unit mass in cgs (U_M)"]
unit_time_in_cgs = sim["/Units"].attrs["Unit time in cgs (U_t)"]
unit_temp_in_cgs = sim["/Units"].attrs["Unit temperature in cgs (U_T)"]
unit_vel_in_cgs = unit_length_in_cgs / unit_time_in_cgs
unit_energy_in_cgs = unit_mass_in_cgs * unit_vel_in_cgs * unit_vel_in_cgs
unit_length_in_si = 0.01 * unit_length_in_cgs
unit_mass_in_si = 0.001 * unit_mass_in_cgs
unit_time_in_si = unit_time_in_cgs
unit_density_in_cgs = unit_mass_in_cgs*unit_length_in_cgs**-3
unit_pressure_in_cgs = unit_mass_in_cgs/unit_length_in_cgs*unit_time_in_cgs**-2
unit_int_energy_in_cgs = unit_energy_in_cgs/unit_mass_in_cgs
unit_entropy_in_cgs = unit_energy_in_cgs/unit_temp_in_cgs
Myr_in_cgs = 3.154e13

pos = sim["/PartType0/Coordinates"][:,:]
x = pos[:,0] - boxSize / 2
y = pos[:,1] - boxSize / 2
z = pos[:,2] - boxSize / 2
vel = sim["/PartType0/Velocities"][:,:]
r = sqrt(x**2 + y**2 + z**2)
v_r = (x * vel[:,0] + y * vel[:,1] + z * vel[:,2]) / r
u = sim["/PartType0/InternalEnergy"][:]
S = sim["/PartType0/Entropy"][:]
P = sim["/PartType0/Pressure"][:]
rho = sim["/PartType0/Density"][:]

# Bin te data
r_bin_edge = np.arange(0., 0.5, 0.01)
r_bin = 0.5*(r_bin_edge[1:] + r_bin_edge[:-1])
rho_bin,_,_ = stats.binned_statistic(r, rho, statistic='mean', bins=r_bin_edge)
v_bin,_,_ = stats.binned_statistic(r, v_r, statistic='mean', bins=r_bin_edge)
P_bin,_,_ = stats.binned_statistic(r, P, statistic='mean', bins=r_bin_edge)
S_bin,_,_ = stats.binned_statistic(r, S, statistic='mean', bins=r_bin_edge)
u_bin,_,_ = stats.binned_statistic(r, u, statistic='mean', bins=r_bin_edge)
rho2_bin,_,_ = stats.binned_statistic(r, rho**2, statistic='mean', bins=r_bin_edge)
v2_bin,_,_ = stats.binned_statistic(r, v_r**2, statistic='mean', bins=r_bin_edge)
P2_bin,_,_ = stats.binned_statistic(r, P**2, statistic='mean', bins=r_bin_edge)
S2_bin,_,_ = stats.binned_statistic(r, S**2, statistic='mean', bins=r_bin_edge)
u2_bin,_,_ = stats.binned_statistic(r, u**2, statistic='mean', bins=r_bin_edge)
rho_sigma_bin = np.sqrt(rho2_bin - rho_bin**2)
v_sigma_bin = np.sqrt(v2_bin - v_bin**2)
P_sigma_bin = np.sqrt(P2_bin - P_bin**2)
S_sigma_bin = np.sqrt(S2_bin - S_bin**2)
u_sigma_bin = np.sqrt(u2_bin - u_bin**2)


# Plot the interesting quantities
figure()

# Velocity profile --------------------------------
subplot(231)
plot(r * unit_length_in_cgs, v_r * unit_vel_in_cgs, '.', color='r', ms=0.5, alpha=0.4)
#errorbar(r_bin, v_bin, yerr=v_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Radius}}~r (cm)$", labelpad=0)
ylabel("${\\rm{Radial~velocity}}~v_r (cm \cdot s^{-1})$", labelpad=0)
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
xlim(0, 3e22)
ylim(-2e7, 4e7)

# Density profile --------------------------------
subplot(232)
plot(r * unit_length_in_cgs, rho * unit_density_in_cgs, '.', color='r', ms=0.5, alpha=0.4)
#errorbar(r_bin, rho_bin, yerr=rho_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Radius}}~r (cm)$", labelpad=0)
ylabel("${\\rm{Density}}~\\rho (g \cdot cm^{-3})$", labelpad=2)
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
xlim(0, 3e22)
ylim(0, 5e24)

# Pressure profile --------------------------------
subplot(233)
plot(r * unit_length_in_cgs, P * unit_pressure_in_cgs, '.', color='r', ms=0.5, alpha=0.4)
#errorbar(r_bin, P_bin, yerr=P_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Radius}}~r (cm)$", labelpad=0)
ylabel("${\\rm{Pressure}}~P (Ba)$", labelpad=0)
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
xlim(0, 3e22)
ylim(0, 1e-10)

# Internal energy profile -------------------------
subplot(234)
plot(r * unit_length_in_cgs, u * unit_int_energy_in_cgs, '.', color='r', ms=0.5, alpha=0.4)
#errorbar(r_bin, u_bin, yerr=u_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Radius}}~r (cm)$", labelpad=0)
ylabel("${\\rm{Internal~Energy}}~u (erg \cdot g^{-1})$", labelpad=0)
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
xlim(0, 3e22)
ylim(0, 1e16)

# Entropy profile ---------------------------------
subplot(235)
plot(r * unit_length_in_cgs, S * unit_energy_in_cgs, '.', color='r', ms=0.5, alpha=0.4)
#errorbar(r_bin, S_bin, yerr=S_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Radius}}~r (cm)$", labelpad=0)
ylabel("${\\rm{Entropy}}~S (erg \cdot K^{-1})$", labelpad=0)
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
xlim(0, 3e22)
ylim(0, 1e56)

# Information -------------------------------------
subplot(236, frameon=False)

text(-0.49, 0.7, "Time (Myr) $%.2f$"%(time*unit_time_in_cgs/Myr_in_cgs), fontsize=10)
xlim(-0.5, 0.5)
ylim(0, 1)
xticks([])
yticks([])


savefig("sedov.png", dpi=200)





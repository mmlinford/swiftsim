"""
Plots the mean temperature.
"""

import matplotlib.pyplot as plt

from swiftsimio import load

from unyt import Myr, K, mh, cm

import numpy as np


def setup_axes():
    """
    Sets up the axes and returns fig, ax
    """

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(4, 8), sharex=True)

    ax = ax.flatten()

    ax[0].semilogy()

    ax[1].set_xlabel("Simulation time elapsed [Myr]")
    ax[0].set_ylabel("Temperature of Universe [K]")
    ax[1].set_ylabel("Physical Density of Universe [$n_H$ cm$^{-3}$]")

    ax[0].set_xlim(0, 50)

    fig.tight_layout()

    return fig, ax


def get_data(handle: float, n_snaps: int):
    """
    Returns the elapsed simulation time, temperature, and density
    """

    t0 = 0.0

    times = []
    temps = []
    densities = []

    for snap in range(n_snaps):
        data = load(f"{handle}_{snap:04d}.hdf5")

        if snap == 0:
            t0 = data.metadata.t.to(Myr).value

        times.append(data.metadata.t.to(Myr).value - t0)
        temps.append(np.mean(data.gas.temperature.to(K).value))
        densities.append(np.mean(data.gas.density.to(mh / cm**3).value) / (data.metadata.scale_factor**3) )

    return times, temps, densities


def plot_single_data(handle: str, n_snaps: int, label: str, ax: plt.axes):
    """
    Takes the a single simulation and plots it on the axes.
    """

    data = get_data(handle, n_snaps)

    ax[0].plot(
        data[0], data[1],
        label=label,
        marker="o",
        ms=2
    )

    ax[1].plot(
        data[0], data[2],
        label=label,
        marker="o",
        ms=2
    )

    
    return


def make_plot(handles, names, n_snaps=100):
    """
    Makes the whole plot and returns the fig, ax objects.
    """

    fig, ax = setup_axes()

    for handle, name in zip(handles, names):
        plot_single_data(handle, n_snaps, name, ax)

    ax[0].legend()

    return fig, ax 


if __name__ == "__main__":
    """
    Plot everything!
    """

    handles = ["redshift_dependence_no_z", "redshift_dependence_low_z"]
    names = ["No Cosmology", "Low Redshift"]

    fig, ax = make_plot(handles, names)

    fig.savefig("redshift_dependence.png", dpi=300)


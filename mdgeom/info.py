# -*- coding: utf-8 -*-
"""
Information about simulations
=============================

This module contains functions that provide overview information about a
trajectory that may be useful for a Methods section of a paper.

"""
import prettytable


def summary(*universes, labels=None):
    """Print summary information about all `universes`.

    Parameters
    ----------
    universes : list of Universe
        Print information about these simulation systems
        in a table.

    labels: list of str, optional
        Add a label to the output, e.g., the name or ID of the simulation.
        There needs to be one element on `labels` for each universe.

    """
    data_column_keys = [
        "n_atoms",
        "Lx",
        "Ly",
        "Lz",
        "alpha",
        "beta",
        "gamma",
        "n_frames",
        "totaltime",
        "dt",
    ]
    data_column_names = [
        "n_atoms",
        "Lx/Å",
        "Ly/Å",
        "Lz/Å",
        "alpha",
        "beta",
        "gamma",
        "n_frames",
        "totaltime/ns",
        "dt/ps",
    ]

    if len(labels) != len(universes):
        raise ValueError(
            f"labels {labels} must contain one "
            f"label per universe {universes}"
        )

    if labels is None:
        labels = len(universes) * [None]
        column_names = data_column_names
    else:
        column_names = ["simulation"] + data_column_names

    table = prettytable.PrettyTable()
    table.field_names = column_names

    for u, label in zip(universes, labels):
        data = extract(u)
        data["totaltime"] /= 1000  # convert ps to ns
        row = [data[key] for key in data_column_keys]
        if label is not None:
            row = [label] + row
        table.add_row(row)

    print(table)


def extract(u):
    """Extract summary information from a universe.

    The information is returned as a dictionary for keys
    - `n_atoms`: number of atoms
    - `Lx`, `Ly`, `Lz`: length of the unit cell in Å (from first frame of the
      trajectory)
    - `alpha`, `beta`, `gamma`: angles of the unit cell in degrees
    - `n_frames`: number of frames in the trajectory
    - `totaltime`: total simulation time in ps
    - `dt`: time between saved frames in the trajectory in ps


    Parameters
    ----------
    u : Universe
        MDAnalysis Universe.

    Returns
    -------
    data : dict

    """

    try:
        Lx, Ly, Lz, alpha, beta, gamma = u.dimensions
    except TypeError:
        # universe without a regular box
        Lx, Ly, Lz, alpha, beta, gamma = 0, 0, 0, 0, 0, 0
    data = {
        "n_atoms": u.atoms.n_atoms,
        "Lx": Lx,
        "Ly": Ly,
        "Lz": Lz,
        "alpha": alpha,
        "beta": beta,
        "gamma": gamma,
        "n_frames": u.trajectory.n_frames,
        "totaltime": u.trajectory.totaltime,
        "dt": u.trajectory.dt,
    }
    return data

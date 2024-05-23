import numpy as np
from numpy.testing import assert_allclose
import pytest

from mdgeom import distances

# ADK FRET pairs
# see doi: 10.1016/j.jmb.2009.09.009.
fret_pairs = [
    ("resname ALA and resid 55", "resname VAL and resid 169"),
    ("resname ALA and resid 127", "resname ALA and resid 194"),
    ("resname ILE and resid 52", "resname LYS* and resid 145"),
]


@pytest.fixture
def u_nopbc():
    # simulation without periodic boundaries
    # (box == None)
    return mda.Universe(data.PSF, data.DCD)


@pytest.fixture
def u_pbc():
    # simulation with periodic boundaries and protein split
    # across unit cell
    return mda.Universe(data.TPR, data.XTC)


#
# Tests for distances.minimage_distance()
#


def _minimage_distance(u):
    # helper function
    a = sum(u.select_atoms(pair[0] + " and name CA") for pair in fret_pairs)
    b = sum(u.select_atoms(pair[1] + " and name CA") for pair in fret_pairs)

    r = distances.minimage_distance(a, b)

    return r


def test_minimage_distance_nopbc(u_nopbc):
    reference = np.array([12.91501288, 30.57278011, 28.92306815])
    r = _minimage_distance(u_nopbc)
    assert_allclose(r, reference)


def test_minimage_distance_pbc(u_pbc):
    reference = np.array([29.67992244, 38.36642607, 42.18877151])
    r = _minimage_difference_vector(u_pbc)
    assert_allclose(r, reference)


#
# Tests for distances.minimage_difference_vector()
#


def _minimage_difference_vector(u):
    # helper function
    a = sum(u.select_atoms(pair[0] + " and name CA") for pair in fret_pairs)
    b = sum(u.select_atoms(pair[1] + " and name CA") for pair in fret_pairs)

    r = distances.minimage_difference_vector(a, b)

    return r


def test_minimage_difference_vector_nopbc(u_nopbc):
    reference = np.array(
        [
            [-6.7521524, -7.5131125, -8.047306],
            [-21.687479, -15.805498, -14.646992],
            [16.235527, -19.996527, 13.156384],
        ]
    )
    r = _minimage_difference_vector(u_nopbc)
    assert_allclose(r, reference)


def test_minimage_difference_vector_pbc(u_pbc):
    reference = np.array(
        [
            [-11.659996, 8.610001, 25.900003],
            [-9.439999, 37.177002, 0.8600006],
            [-31.1585, 21.841507, -18.220566],
        ]
    )
    r = _minimage_difference_vector(u_pbc)
    assert_allclose(r, reference)

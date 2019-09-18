"""Utils.py."""

from __future__ import division, print_function

import astropy.constants as const
import astropy.units as u
import scipy as sp


def energy(l):
    """Calculate the energy of a photon with wavelength l.

    Parameters
    ----------
    l : float
        Photon's wavelength in meters.

    Returns
    -------
    E : float
        Energy.

    """
    h = const.h
    c = const.c
    E = h * c / l
    return E


def airmass(z):
    """Calculate airmass given zenith angle.

    Reference:
    R.H. Hardie, 1962, 'Photoelectric Reductions', Chapter 8 of Astronomical
    Techniques, W.A. Hiltner (Ed), Stars and Stellar Systems, II
    (University of Chicago Press: Chicago), pp178-208.

    Parameters
    ----------
    z : float
        Zenith angle.

    Returns
    -------
    X : float
        Airmass.

    """
    DEG_TO_RAD = sp.pi / 180
    secz = 1 / sp.cos(z * DEG_TO_RAD)
    X = secz - 0.0018167 * (secz - 1) - 0.002875 * \
        (secz - 1) ** 2 - 0.0008083 * (secz - 1) ** 3
    return X


def flux(m, Z):
    """Calculate total flux from a source.

    Uses the formula: m = -2.5 * log(f) + Z

    Parameters
    ----------
    m : float
        Airmass, extinction corrected apparent magnitude.
    Z : float
        Zero-point constant for flux.

    Returns
    -------
    f : float
        The flux from the object.

    """
    f = 10 ** (-.4 * m) * Z
    return f * u.erg / u.s / u.cm**2 / u.angstrom


def correct_mag(m0, X, k):
    """Correct magnitude for airmass and extinction.

    Parameters
    ----------
    m0 : float
        Apparent magnitude.
    X : float
        Airmass.
    k : float
        Extinction coefficient.

    Returns
    -------
    m : float
        The corrected apparent magnitude.

    """
    m = m0 + k * X
    return m


def s_n(n, N, phi, m, dark, t):
    """Calculate Signal-to-noise.

    Given the number of photons from the background and the source, the
    number of pixels  and the quantum efficiency it calculates the
    signal-to-noise.

    Parameters
    ----------
    n : float
        Photons from source.
    N : float
        Photons from background.
    phi : float
        Read-out noise.
    m : int
        Number of pixels.
    dark : float
        Dark current.
    t : float
        Integration time

    Returns
    -------
    s_n : float
        The signal-to-noise.

    """
    sn = n * t / sp.sqrt(n * t + N * t * m + m * dark * t + m * phi**2)
    return sns

"""Wave shoaling adjustment for significant wave height."""

from typing import Union

import numpy as np

from cht_utils.physics.disper import disper


def deshoal(
    hm0: Union[float, np.ndarray],
    Tp: float,
    d_profile: Union[float, np.ndarray],
    d_BC: float,
) -> Union[float, np.ndarray]:
    """Adjust significant wave height for shoaling.

    Parameters
    ----------
    hm0 : float or np.ndarray
        Significant wave height at the boundary.
    Tp : float
        Peak wave period.
    d_profile : float or np.ndarray
        Water depth(s) along the profile.
    d_BC : float
        Water depth at the boundary condition location.

    Returns
    -------
    float or np.ndarray
        De-shoaled significant wave height.
    """
    cg_profile = wavecelerity(Tp, d_profile)
    cg_BC = wavecelerity(Tp, d_BC)
    return hm0 * np.sqrt(cg_profile / cg_BC)


def wavecelerity(
    Tp: float,
    d: Union[float, np.ndarray],
    g: float = 9.81,
) -> Union[float, np.ndarray]:
    """Compute wave group velocity.

    Parameters
    ----------
    Tp : float
        Peak wave period.
    d : float or np.ndarray
        Water depth.
    g : float
        Gravitational acceleration.

    Returns
    -------
    float or np.ndarray
        Group velocity.
    """
    k = disper(2 * np.pi / Tp, d, g)
    n = 0.5 * (1 + 2 * k * d / np.sinh(2 * k * d))
    c = g * Tp / (2 * np.pi) * np.tanh(k * d)
    return n * c

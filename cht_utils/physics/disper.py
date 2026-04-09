"""Linear wave dispersion relation solver.

Absolute error in k*h < 5.0e-16 for all k*h.

Example::

    k = disper(2 * np.pi / 5, 5, 9.81)
"""

from typing import List, Union

import numpy as np


def disper(
    w: Union[float, List[float]],
    h: Union[float, List[float]],
    g: float = 9.81,
) -> np.ndarray:
    """Solve the linear dispersion relation for wave number.

    Parameters
    ----------
    w : float or list of float
        Angular frequency (2*pi/T).
    h : float or list of float
        Water depth.
    g : float
        Gravitational acceleration.

    Returns
    -------
    np.ndarray
        Wave number(s).
    """
    if not isinstance(w, list):
        w = [w]
    if not isinstance(h, list):
        h = [h]

    w2 = [iw**2 * ih / g for iw, ih in zip(w, h)]
    q = [iw2 / (1 - np.exp(-(iw2 ** (5 / 4)))) ** (2 / 5) for iw2 in w2]

    the = np.tanh(q)
    the2 = 1 - the**2

    a = (1 - q * the) * the2
    b = the + q * the2
    c = q * the - w2

    D = b**2 - 4 * a * c
    arg = (-b + np.sqrt(D)) / (2 * a)
    iq = np.where(D < 0)[0]
    if len(iq) > 0:
        arg[iq] = -c[iq] / b[iq]
    q = q + arg

    k = np.sign(w) * q / h
    if np.isnan(k).any():
        k = np.array(k)
        k[np.isnan(k)] = 0

    return k


def disper_fentonmckee(
    sigma: Union[float, np.ndarray],
    d: Union[float, np.ndarray],
    g: float = 9.81,
) -> Union[float, np.ndarray]:
    """Fenton-McKee approximation of the linear dispersion relation.

    Parameters
    ----------
    sigma : float or np.ndarray
        Angular frequency.
    d : float or np.ndarray
        Water depth.
    g : float
        Gravitational acceleration.

    Returns
    -------
    float or np.ndarray
        Approximate wave number(s).
    """

    def coth(x):
        return 1.0 / np.tanh(x)

    return (sigma**2 / g) * (coth(sigma * np.sqrt(d / g)) ** (3 / 2)) ** (2 / 3)

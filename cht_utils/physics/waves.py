"""Wave decomposition utilities."""

import numpy as np
import pandas as pd
from scipy.signal import detrend


def split_waves_guza(df: pd.DataFrame, zb: float) -> pd.DataFrame:
    """Decompose water surface into incident and reflected waves (Guza method).

    Performs a time-domain decomposition of the water surface elevation and
    velocity into incoming and outgoing components using shallow-water theory.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns ``"z"`` (water level) and ``"u"`` (velocity).
    zb : float
        Bed level elevation.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with added columns ``"zin"``, ``"zout"``,
        ``"uin"``, ``"uout"``.
    """
    g = 9.81

    zsi = df["z"].values
    umi = df["u"].values

    zsm = np.mean(zsi)
    umm = np.mean(umi)
    h = zsm - zb

    zs = zsi - zsm
    um = umi - umm

    zsd = detrend(zs, type="linear")
    umd = detrend(um, type="linear")

    zsm = zsm + (zs - zsd)
    zs = zsd
    umm = umm + (um - umd)
    um = umd

    hh = zs + h
    c = np.sqrt(g * hh)
    q = umd * hh

    ein = (zs * c + q) / (2 * c)
    eout = (zs * c - q) / (2 * c)

    df["zin"] = ein + zsm
    df["zout"] = eout + zsm
    df["uin"] = np.sqrt(1.0 / hh**2) * c * ein + umm
    df["uout"] = -np.sqrt(1.0 / hh**2) * c * eout + umm

    return df

'''
DISPER  Linear dispersion relation.

absolute error in k*h < 5.0e-16 for all k*h

Syntax:
k = disper(omega, h, g)

Input:
w = 2*pi/T, where T is the wave period
h = water depth
g = gravity constant

Output:
k = wave number

Example
k = disper(2*pi/5,5,9.81);
'''

import numpy as np
from scipy.signal import detrend

def disper(w, h, g=9.81):
    if not type(w) == list:
        w=[w]
    if not type(h) == list:
        h=[h]
    w2 = [iw**2 * ih / g for (iw,ih) in zip(w,h)]
    q = [iw2 / (1 - np.exp(-(iw2**(5/4))))**(2/5) for iw2 in w2]

    thq = np.tanh(q)
    thq2 = 1 - thq**2

    a = (1 - q * thq) * thq2
    b = thq + q * thq2
    c = q * thq - w2

    D = b**2 - 4 * a * c
    arg = (-b + np.sqrt(D)) / (2 * a)
    iq = np.where(D < 0)[0]
    if iq:
        print(iq)
        arg[iq] = -c[iq] / b[iq]
    q = q + arg

    k = np.sign(w) * q / h
    if np.isnan(k).any():
        k = np.array(k)
        k[np.isnan(k)] = 0

    return k

def disper_fentonmckee(sigma, d, g=9.81):

    def coth(x):
        return 1. / np.tanh(x)

    k = (sigma**2 / g) * (coth(sigma * np.sqrt(d / g))**(3 / 2))**(2 / 3)
    # k = sigma**2 / g * coth((sigma * np.sqrt(d / g))**(3 / 2))**(2 / 3)

    return k

def split_waves_guza(df, zb):
    """Guza split waves"""

    #  in timespace

    g = 9.81

    zsi = df["z"].values
    umi = df["u"].values

    # Mean water level and velocity
    zsm = np.mean(zsi)
    umm = np.mean(umi)

    # Mean water depth
    h = zsm - zb

    # Adjust to zero-centered water level and velocity
    zs = zsi - zsm
    um = umi - umm

    zsd = detrend(zs, type="linear")
    umd = detrend(um, type="linear")

    zsm = zsm + (zs - zsd)
    zs  = zsd
    umm = umm + (um - umd)
    um = umd

    hh = zs + h
    c = np.sqrt(g * hh)
    q = umd * hh
    ein = (zs * c + q) / (2 * c)
    eout = (zs * c - q) / (2 * c)
    zsin = ein + zsm
    zsout = eout + zsm
    uin   = (np.sqrt(1.0 / (hh**2)) * c * ein) + umm
    uout  = - (np.sqrt(1.0/(hh**2)) * c * eout) + umm

    df["zin"]  = zsin
    df["zout"] = zsout
    df["uin"]   = uin
    df["uout"]  = uout

    return df

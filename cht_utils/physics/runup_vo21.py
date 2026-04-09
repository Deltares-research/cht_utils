"""Wave runup prediction using Van Ormondt (2021) empirical formulation.

Computes the 2% exceedance runup level from offshore wave conditions
and nearshore profile information, decomposed into setup, low-frequency,
and high-frequency components with directional spreading corrections.
"""

from typing import Any, Union

import numpy as np
from scipy import interpolate as intp

from cht_utils.physics.disper import disper, disper_fentonmckee


class runup_vo21:
    """Wave runup calculator following Van Ormondt (2021).

    Parameters
    ----------
    hm0 : np.ndarray
        Significant wave height.
    tp : float
        Peak wave period.
    ztide : float or np.ndarray
        Tidal water level.
    sl1 : float, np.ndarray, or dict
        Surf slope parameter. Interpretation depends on *sl1opt*.
    sl2 : float or np.ndarray
        Beach slope (foreshore).
    sl1opt : str
        Slope option: ``"ad"`` (adaptive), ``"slope"`` (direct), or
        ``"xz"`` (cross-shore profile dict with keys ``"x"`` and ``"z"``).
    """

    def __init__(
        self,
        hm0: np.ndarray,
        tp: float,
        ztide: Union[float, np.ndarray],
        sl1: Any,
        sl2: Union[float, np.ndarray],
        sl1opt: str,
    ) -> None:
        self.hm0 = hm0
        self.tp = tp
        self.h = self.tp * 10

        k = disper(2 * np.pi / self.tp, self.h, 9.81)
        self.l1 = 2 * np.pi / k
        self.steepness = self.hm0 / self.l1

        if np.size(ztide) == 1:
            ztide = np.zeros(np.size(self.hm0)) + ztide
        if np.size(sl2) == 1:
            sl2 = np.zeros(np.size(self.hm0)) + sl2

        if sl1opt == "ad":
            gambr = 1
            fh = 1 / gambr
            surfslope2 = np.zeros(np.size(ztide))
            surfslope2[0] = (fh * self.hm0) / (fh * self.hm0 / self.surfslope) ** 1.5
            for itide in range(np.size(ztide)):
                if ztide[itide] == 0:
                    pass
                else:
                    xxx = np.arange(-1000, 5005, 5)
                    yyy0 = -sl1[itide] * xxx ** (2 / 3)
                    yyy0[np.argwhere(xxx < 0)] = -xxx[np.argwhere(xxx < 0)] * sl2[itide]
                    yyy = yyy0 - ztide[itide]
                    ii1 = np.argwhere(yyy < 0)[0]
                    ii2 = np.argwhere(yyy < -self.hm0[itide] / gambr)[0]
                    surfslope2[itide] = (yyy[ii1] - yyy[ii2]) / (xxx[ii2] - xxx[ii1])
        elif sl1opt == "slope":
            surfslope2 = sl1
        elif sl1opt == "xz":
            gambr = 1
            xxxx = np.arange(sl1["x"][0], sl1["x"][-1], 1)
            zzzz = intp.interp1d(sl1["x"], sl1["z"])(xxxx)
            sl1["x"] = xxxx
            sl1["z"] = zzzz
            surfslope2 = np.zeros(np.size(ztide))
            xxx = sl1["x"]
            yyy0 = sl1["z"]
            for itide in range(np.size(ztide)):
                yyy = yyy0 - ztide[itide]
                ii1 = np.argwhere(yyy < 0)[0]
                ii2 = np.argwhere(yyy < -self.hm0[itide] / gambr)[0]
                surfslope2[itide] = (yyy[ii1] - yyy[ii2]) / (xxx[ii2] - xxx[ii1])

        self.sl2 = sl2
        self.ztide = ztide
        self.surfslope = surfslope2
        self.ksis = self.surfslope / np.sqrt(self.hm0 / self.l1)

    def compute_r2p(self, drspr: np.ndarray, phi: np.ndarray) -> None:
        """Compute the 2% exceedance runup level.

        Parameters
        ----------
        drspr : np.ndarray
            Directional spreading (degrees).
        phi : np.ndarray
            Wave direction relative to shore normal (degrees).
        """
        ksi1 = self.surfslope / np.sqrt(self.hm0 / self.l1)
        ksi2 = self.sl2 / np.sqrt(self.hm0 / self.l1)

        self.compute_setup(drspr, phi)
        self.compute_hm0_lf(drspr, phi)
        self.compute_hm0_hf(drspr, phi)

        self.r2p = self.setup + 0.82396 * np.sqrt(
            (0.82694 * self.hm0_lf) ** 2 + (0.73965 * self.hm0_hf) ** 2
        ) * ksi2**0.15201 * ksi1**(-0.086635)

    def compute_setup(self, drspr: np.ndarray, phi: np.ndarray) -> None:
        """Compute wave setup component."""
        b = [4.0455506e00, 2.3740615e-02, 2.0340287e00, 4.6497588e-01,
             7.0244541e-01, 5.0000000e-01, -3.1727583e-01]
        ksib = self.sl2 / np.sqrt(self.hm0 / self.l1)
        v = b[0] * self.hm0 * (
            b[1] + np.exp(-b[2] * self.ksis ** b[6] * ksib ** b[3]) * ksib ** b[4]
        )
        fac1 = self._dirspreadfac_setup(drspr)
        fac2 = self._directionfac_setup(phi)
        self.setup = v * fac1 * fac2

    def compute_hm0_lf(self, drspr: np.ndarray, phi: np.ndarray) -> None:
        """Compute low-frequency wave height component."""
        b = [3.4547125e00, 5.8790748e-01, 3.6906975e00, 2.3378556e-01,
             2.3038164e00, 0, 5.0000000e-01, 0]

        self.tm0_ig = self._compute_tm01_ig()
        self.IG = self._compute_ig_in()

        l0 = np.squeeze(np.sqrt(9.81 * 0.33333 * self.hm0) * self.tm0_ig)
        ksib = np.squeeze(self.sl2 / np.sqrt(self.IG / l0))
        ksib = np.maximum(ksib - b[5], 0)
        ksibm = b[4] * self.surfslope ** b[3]
        psibd = b[0] * ksib ** b[1]
        psibr = (b[0] * ksibm ** b[1] - 2.0) * np.exp(-(ksib - ksibm) / b[2]) + 2
        psib = psibd.copy()
        ind = np.argwhere(ksib > ksibm)
        psib[ind] = psibr[ind]

        v = self.IG * psib
        fac1 = self._dirspreadfac_lf(drspr)
        fac2 = self._directionfac_lf(phi)
        self.hm0_lf = v * fac1 * fac2

    def compute_hm0_hf(self, drspr: np.ndarray, phi: np.ndarray) -> None:
        """Compute high-frequency wave height component."""
        b = [9.5635099e-01, 2.0143005e00, 5.3602429e-01, 2.0000000e00,
             0.0000000e00, 6.1856544e-01, 1.0000000e00, 0.0000000e00]
        ksib = self.sl2 / np.sqrt(self.hm0 / self.l1)
        self.hm0_hf = (
            b[0] * self.hm0 * ksib ** b[1]
            * np.tanh((self.ksis + b[4]) ** b[5] / (b[2] * ksib ** b[3]))
        )

    def _compute_tm01_ig(self) -> np.ndarray:
        """Compute infragravity wave period."""
        b = [4.4021341e-07, 1.8635421e00, -4.2705433e-01, 7.2541023e-02, 2.0058478e01]
        tm0_ig = (
            b[0] + b[1] * self.surfslope ** b[2] * self.steepness ** b[3]
            + b[4] * self.surfslope
        )
        return tm0_ig * self.tp

    def _compute_ig_in(self) -> np.ndarray:
        """Compute incoming infragravity wave amplitude."""
        b = [2.2740842e00, 1.0, 0.5, 2.7211454e03, 2.0, 1.7794945e01, 1.8728433e-01]
        return self.hm0 * (
            b[0] * self.surfslope ** b[2] * np.exp(-b[5] * self.surfslope ** b[1])
            + b[6] * np.exp(-b[3] * self.steepness ** b[4])
        )

    def _directionfac_setup(self, phi: np.ndarray) -> np.ndarray:
        b = [1.4291, 0.0035124, 0.31891]
        return 1 - b[1] * self.steepness ** b[2] * np.absolute(phi) ** b[0]

    def _dirspreadfac_setup(self, drspr: np.ndarray) -> np.ndarray:
        b = [0.031448, 0.69432, 0.66677]
        return np.exp(-b[0] * drspr ** b[1] * (1.0 - np.tanh(b[2] * self.ksis)))

    def _dirspreadfac_lf(self, drspr: np.ndarray) -> np.ndarray:
        b = [0.047593, 0.67228, 0.50777]
        return np.exp(-b[0] * drspr ** b[1] * (1.0 - np.tanh(b[2] * self.ksis)))

    def _directionfac_lf(self, phi: np.ndarray) -> np.ndarray:
        b = [0.40488, 2.7073]
        return np.cos(phi * np.pi / 180) ** (b[0] + b[1] * self.surfslope)

    def _directionfac_hf(self, phi: np.ndarray) -> np.ndarray:
        b = [4.1355e-10]
        return np.cos(phi * np.pi / 180) ** b[0]

    def _dirspreadfac_hf(self, drspr: np.ndarray) -> np.ndarray:
        b = [0.044544, 1, 21.1281]
        return np.exp(-b[0] * drspr ** b[1] * (1.0 - np.tanh(b[2] * self.ksis)))

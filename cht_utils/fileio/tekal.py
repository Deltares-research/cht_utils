"""Read and write Delft3D Tekal ASCII files.

A Tekal file consists of consecutive blocks, each containing:

- Optional ``*``-prefixed comment lines
- A block name line
- A size line (nrows, ncols, and optionally nx, ny)
- Data lines

Example usage::

    import cht_utils.fileio.tekal as tek

    D = tek.tekal("myfile.tek")
    D.info()
    print(D)
    M = D.read(0)  # load first block
"""

from typing import IO, List, Optional

import numpy as np


class tekalblock:
    """A single data block within a Tekal file."""

    def __init__(self) -> None:
        self.header: List[str] = []
        self.nhdr: Optional[int] = None
        self.index: Optional[int] = None
        self.name: Optional[str] = None
        self.npts: Optional[int] = None
        self.nvar: Optional[int] = None
        self.nx: Optional[int] = None
        self.ny: Optional[int] = None
        self.data: Optional[np.ndarray] = None
        self.tell: float = 0.0

    def __str__(self) -> str:
        return (
            f"{self.index:5g}: nvar={self.nvar:4g} npts={self.npts:5g} "
            f"({self.nx:5g} x{self.ny:4g}) {self.nhdr:4} '{self.name}'\n"
        )

    def strheader(self) -> str:
        """Return a header line describing the column layout."""
        return (
            f"{'block':>5s}:      {'nvar':>4s}      {'npts':>5s} "
            f"({'nx':>5s} x{'ny':>4s}) {'nhdr':>4s} '{'blockname'}'\n"
        )

    def info(self, fid: IO, index: int) -> Optional["tekalblock"]:
        """Read block metadata (header, name, dimensions) from an open file.

        Parameters
        ----------
        fid : IO
            Open file handle.
        index : int
            Block index.

        Returns
        -------
        tekalblock or None
            Self if successful, ``None`` at end-of-file.
        """
        self.tell = fid.tell()

        rec = fid.readline()
        if len(rec) == 0:
            return None

        while rec[0] == "*":
            self.header.append(rec)
            rec = fid.readline()

        self.index = index
        self.name = rec.split()[0]
        self.header.append(rec)
        self.nhdr = len(self.header)

        rec = fid.readline()
        parts = rec.split()
        self.npts = int(parts[0])
        self.nvar = int(parts[1])
        if len(parts) == 2:
            self.nx = self.npts
            self.ny = 1
        else:
            self.nx = int(parts[2])
            self.ny = int(parts[3])

        # Skip data lines
        for _ in range(self.npts):
            fid.readline()

        return self

    def load(self, fid: IO) -> np.ndarray:
        """Load the data from this block into a 3-D numpy array.

        The returned array has shape ``(nvar, nx, ny)``.

        Parameters
        ----------
        fid : IO
            Open file handle, seeked to ``self.tell``.

        Returns
        -------
        np.ndarray
            Data array of shape ``(nvar, nx, ny)``.
        """
        # Skip header
        rec = fid.readline()
        while rec[0] == "*":
            rec = fid.readline()

        rec = fid.readline()
        parts = rec.split()
        npts = int(parts[0])
        nvar = int(parts[1])
        if len(parts) == 2:
            nx, ny = npts, 1
        else:
            nx, ny = int(parts[2]), int(parts[3])

        M = np.zeros([nvar, nx, ny])
        for ix in range(nx):
            for iy in range(ny):
                rec = fid.readline()
                values = rec.split()
                M[:, ix, iy] = [float(v) for v in values[:nvar]]

        self.data = M
        return M

    def check(self) -> None:
        """Assert that all required fields are populated."""
        assert self.header != []
        assert self.nhdr is not None
        assert self.name is not None
        assert self.npts is not None
        assert self.nvar is not None
        assert self.nx is not None
        assert self.ny is not None
        assert self.data is not None

    def copy(self) -> "tekalblock":
        """Create a deep copy of this block."""
        cp = tekalblock()
        cp.header = self.header[:]
        cp.nhdr = self.nhdr
        cp.name = self.name
        cp.npts = self.npts
        cp.nvar = self.nvar
        cp.nx = self.nx
        cp.ny = self.ny
        cp.data = self.data.copy()
        return cp


class tekal:
    """Container for a Tekal file with multiple blocks.

    Parameters
    ----------
    filename : str
        Path to the Tekal file.
    """

    def __init__(self, filename: str) -> None:
        self.filename = filename
        self.blocks: List[tekalblock] = []

    def append(self, block: tekalblock) -> None:
        """Append a validated block to the file structure."""
        block.check()
        self.blocks.append(block)
        self.blocks[-1].index = len(self.blocks) - 1

    def __str__(self) -> str:
        txt = self.blocks[0].strheader()
        for block in self.blocks:
            txt += str(block)
        return txt

    def info(self) -> None:
        """Read metadata for all blocks (without loading data)."""
        fid = open(self.filename, "rb")
        index = 0
        while True:
            B = tekalblock()
            result = B.info(fid, index)
            if result is None or B.npts is None:
                break
            self.blocks.append(B)
            index += 1
        fid.close()

    def read(self, index: int) -> np.ndarray:
        """Load data from a specific block.

        Parameters
        ----------
        index : int
            Block index.

        Returns
        -------
        np.ndarray
            Data array of shape ``(nvar, nx, ny)``.
        """
        fid = open(self.filename, "rb")
        fid.seek(self.blocks[index].tell)
        M = self.blocks[index].load(fid)
        fid.close()
        return M

    def write(self) -> None:
        """Write all blocks to the file."""
        with open(self.filename, "a") as fid:
            for block in self.blocks:
                block.check()
                fid.write("".join(block.header))
                if block.ny == 1:
                    fid.write(f"{block.npts:g} {block.nvar:g}\n")
                else:
                    fid.write(f"{block.npts:g} {block.nvar:g} {block.ny:g}\n")
                for iy in range(block.ny):
                    for ix in range(block.nx):
                        row = " ".join(f"{v:f}" for v in block.data[:, ix, iy])
                        fid.write(row + "\n")

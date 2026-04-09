"""Reader and writer for Deltares-style INI configuration files.

These files consist of ``[Section]`` blocks containing keyword-value pairs
and optional data tables.
"""

import datetime
import re
from typing import Any, List, Optional

import numpy as np
import pandas as pd


class Section:
    """A single section within an INI file.

    Parameters
    ----------
    name : str or None
        Section name (text between square brackets).
    """

    def __init__(self, name: Optional[str] = None) -> None:
        self.name: Optional[str] = name
        self.keyword: List["Keyword"] = []
        self.data: Optional[pd.DataFrame] = None

    def get_value(self, keyword: str) -> Optional[Any]:
        """Return the value of a keyword in this section.

        Parameters
        ----------
        keyword : str
            Keyword name (case-insensitive).

        Returns
        -------
        Any or None
            The keyword value, or ``None`` if not found.
        """
        for kw in self.keyword:
            if kw.name.lower() == keyword.lower():
                return kw.value
        return None


class Keyword:
    """A key-value pair within an INI section.

    Parameters
    ----------
    name : str or None
        Keyword name.
    value : Any
        Keyword value.
    comment : str or None
        Inline comment.
    """

    def __init__(
        self,
        name: Optional[str] = None,
        value: Any = None,
        comment: Optional[str] = None,
    ) -> None:
        self.name = name
        self.value = value
        self.comment = comment


class IniStruct:
    """Deltares INI file parser and writer.

    Parameters
    ----------
    filename : str or None
        If provided, the file is read immediately on construction.
    """

    def __init__(self, filename: Optional[str] = None) -> None:
        self.section: List[Section] = []
        if filename:
            self.read(filename)

    def read(self, filename: str) -> None:
        """Parse an INI file into sections, keywords, and data tables.

        Parameters
        ----------
        filename : str
            Path to the INI file.
        """
        self.section = []
        istart: List[int] = []

        with open(filename, "r") as fid:
            lines = fid.readlines()

        # Find section starts
        for i, line in enumerate(lines):
            ll = line.strip()
            if len(ll) == 0:
                continue
            if ll[0] == "[" and ll[-1] == "]":
                sec = Section(name=ll[1:-1])
                istart.append(i)
                self.section.append(sec)

        # Parse each section
        for isec in range(len(self.section)):
            i1 = istart[isec] + 1
            i2 = istart[isec + 1] - 1 if isec < len(self.section) - 1 else len(lines)

            df = pd.DataFrame()

            for iline in range(i1, i2):
                ll = lines[iline].strip()
                if len(ll) == 0 or ll[0] == "#":
                    continue

                if "=" in ll:
                    key = Keyword()

                    # Handle embedded # in values (e.g. hex colours)
                    if "#" in ll:
                        ipos = [m.start() for m in re.finditer("#", ll)]
                        if len(ipos) > 1:
                            ll = (
                                ll[: ipos[0]]
                                + ll[ipos[0] + 1 : ipos[1]]
                                + ll[ipos[1] + 1 :]
                            )

                    if "#" in ll:
                        j = ll.index("#")
                        key.comment = ll[j + 1 :].strip()
                        ll = ll[:j].strip()

                    tx = ll.split("=")
                    key.name = tx[0].strip()
                    key.value = tx[1].strip()
                    self.section[isec].keyword.append(key)
                else:
                    a_list = ll.split()
                    values = []
                    for item in a_list:
                        try:
                            values.append(float(item))
                        except ValueError:
                            values.append(item)
                    df = pd.concat([df, pd.Series(values)], axis=1)

            if not df.empty:
                df = df.transpose().set_index([0])
                self.section[isec].data = df

    def write(self, file_name: str) -> None:
        """Write the INI structure back to a file.

        Parameters
        ----------
        file_name : str
            Output file path.
        """
        try:
            fid = open(file_name, "w")
        except OSError:
            print(f"Warning! Could not create file: {file_name}")
            return

        try:
            for section in self.section:
                fid.write(f"[{section.name}]\n")

                for kw in section.keyword:
                    valstr = ""
                    if kw.value is not None:
                        if isinstance(kw.value, (float, int)):
                            valstr = str(kw.value)
                        elif isinstance(kw.value, datetime.date):
                            valstr = kw.value.strftime("%Y%m%d")
                        else:
                            valstr = kw.value

                    kwstr = f"   {kw.name:<20s} = "
                    valstr = f"{valstr:<30s}"
                    comstr = f" # {kw.comment}" if kw.comment else ""
                    fid.write(f"{kwstr}{valstr}{comstr}\n")

                if section.data is not None:
                    tt = section.data.index.to_numpy()
                    vv = section.data.array.to_numpy()
                    for it, v in enumerate(vv):
                        if np.isnan(v):
                            v = -999.0
                        fid.write(f"{tt[it]:12.1f}{v:12.3f}\n")

                fid.write("\n")
        finally:
            fid.close()

    def get_value(self, section_name: str, keyword_name: str) -> Optional[Any]:
        """Get a keyword value from a named section.

        Parameters
        ----------
        section_name : str
            Section name.
        keyword_name : str
            Keyword name.

        Returns
        -------
        Any or None
            The keyword value, or ``None`` if not found.
        """
        for section in self.section:
            if section.name == section_name:
                for kw in section.keyword:
                    if kw.name == keyword_name:
                        return kw.value
        return None

    def set_value(
        self,
        section_name: str,
        keyword_name: str,
        keyword_value: Any,
        keyword_comment: Optional[str] = None,
    ) -> None:
        """Set a keyword value, creating the section or keyword if needed.

        Parameters
        ----------
        section_name : str
            Section name.
        keyword_name : str
            Keyword name.
        keyword_value : Any
            New value.
        keyword_comment : str or None
            Optional inline comment.
        """
        for section in self.section:
            if section.name == section_name:
                for kw in section.keyword:
                    if kw.name == keyword_name:
                        kw.value = keyword_value
                        return
                kw = Keyword(keyword_name, keyword_value, keyword_comment)
                section.keyword.append(kw)
                return

        kw = Keyword(keyword_name, keyword_value, keyword_comment)
        section = Section(section_name)
        section.keyword.append(kw)
        self.section.append(section)

    def get_data(
        self,
        section_name: str,
        keyword_list: Optional[List[str]] = None,
        value_list: Optional[List[str]] = None,
    ) -> Optional[pd.DataFrame]:
        """Return the data table from a matching section.

        Parameters
        ----------
        section_name : str
            Section name.
        keyword_list : List[str] or None
            Keywords to match for disambiguation.
        value_list : List[str] or None
            Values that the keywords must have.

        Returns
        -------
        pd.DataFrame or None
        """
        for section in self.section:
            if section.name == section_name:
                if keyword_list and value_list:
                    ok = True
                    for key, val in zip(keyword_list, value_list):
                        for kw in section.keyword:
                            if kw.name == key and kw.value != val:
                                ok = False
                    if not ok:
                        continue
                return section.data
        return None

    def get_section(
        self,
        section_name: str,
        keyword_list: Optional[List[str]] = None,
        value_list: Optional[List[str]] = None,
    ) -> Optional[Section]:
        """Return the first section matching the given name and filters.

        Parameters
        ----------
        section_name : str
            Section name.
        keyword_list : List[str] or None
            Keywords to match for disambiguation.
        value_list : List[str] or None
            Values that the keywords must have.

        Returns
        -------
        Section or None
        """
        for section in self.section:
            if section.name == section_name:
                if keyword_list and value_list:
                    ok = True
                    for key, val in zip(keyword_list, value_list):
                        for kw in section.keyword:
                            if kw.name == key and kw.value != val:
                                ok = False
                    if not ok:
                        continue
                return section
        return None

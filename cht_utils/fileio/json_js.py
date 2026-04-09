"""I/O for JavaScript-wrapped JSON and CSV files.

These formats prepend a JavaScript variable declaration line before the actual
JSON or CSV payload.
"""

import json
import os
from typing import Any, Union


def read_json_js(file_name: str) -> Any:
    """Read a JavaScript-wrapped JSON file (skips the first line).

    Parameters
    ----------
    file_name : str
        Input file path.

    Returns
    -------
    Any
        Parsed JSON object.
    """
    with open(file_name, "r") as f:
        lines = f.readlines()
    jsn_string = "".join(line.strip() for line in lines[1:])
    return json.loads(jsn_string)


def write_json_js(file_name: str, jsn: Union[list, dict], first_line: str) -> None:
    """Write a JSON object to a JavaScript-wrapped file.

    Parameters
    ----------
    file_name : str
        Output file path.
    jsn : list or dict
        JSON-serializable data.
    first_line : str
        JavaScript variable declaration line (e.g. ``"var data ="``)
    """
    with open(file_name, "w") as f:
        f.write(first_line + "\n")
        if isinstance(jsn, list):
            f.write("[\n")
            for ix, item in enumerate(jsn):
                sep = "," if ix < len(jsn) - 1 else ""
                f.write(json.dumps(item) + sep + "\n")
            f.write("]\n")
        else:
            f.write(json.dumps(jsn) + "\n")


def write_csv_js(file_name: str, csv_string: str, first_line: str) -> None:
    """Write a CSV string to a JavaScript-wrapped file.

    Parameters
    ----------
    file_name : str
        Output file path.
    csv_string : str
        CSV content.
    first_line : str
        JavaScript variable declaration line.
    """
    os.makedirs(os.path.dirname(file_name), exist_ok=True)
    csv_string = csv_string.replace(chr(13), "")
    with open(file_name, "w") as f:
        f.write(first_line + "\n")
        f.write(csv_string)
        f.write("`;")

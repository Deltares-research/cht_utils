"""YAML file read/write utilities."""

from typing import Any, Dict, Optional

import yaml


def dict2yaml(file_name: str, dct: Dict[str, Any], sort_keys: bool = False) -> None:
    """Write a dictionary to a YAML file.

    Parameters
    ----------
    file_name : str
        Output file path.
    dct : Dict[str, Any]
        Dictionary to serialize.
    sort_keys : bool
        Whether to sort keys alphabetically.
    """
    yaml_string = yaml.dump(dct, sort_keys=sort_keys)
    with open(file_name, "w") as f:
        f.write(yaml_string)


def yaml2dict(file_name: str) -> Optional[Dict[str, Any]]:
    """Read a YAML file into a dictionary.

    Parameters
    ----------
    file_name : str
        Input file path.

    Returns
    -------
    Dict[str, Any] or None
        Parsed dictionary, or ``None`` for empty files.
    """
    with open(file_name, "r") as f:
        return yaml.load(f, Loader=yaml.FullLoader)

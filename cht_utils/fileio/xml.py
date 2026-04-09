"""XML serialization and deserialization with type conversion.

Supports reading/writing Python objects and dictionaries as XML, with
optional type attributes (``float``, ``int``, ``datetime``) for automatic
value conversion.
"""

import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime
from typing import Any, IO, Optional


class XMLObject:
    """Marker base class for XML-deserialized objects."""

    pass


def write_node(file: IO, node: Any, nspaces: int) -> None:
    """Recursively write a node (dict or object) to an XML file.

    Parameters
    ----------
    file : IO
        Open file handle.
    node : dict or object
        Node to serialize.
    nspaces : int
        Current indentation level.
    """
    node_dict = node if isinstance(node, dict) else node.__dict__
    spaces = " " * nspaces

    for key, children in node_dict.items():
        for child in children:
            if isinstance(child, dict) or (
                not isinstance(child, dict) and not hasattr(child, "value")
            ):
                if hasattr(child, "value"):
                    _write_leaf(file, spaces, key, child)
                else:
                    file.write(f"{spaces}<{key}>\n")
                    write_node(file, child, nspaces + 2)
                    file.write(f"{spaces}</{key}>\n")
            else:
                _write_leaf(file, spaces, key, child)


def _write_leaf(file: IO, spaces: str, key: str, node: Any) -> None:
    """Write a leaf node with optional type attribute."""
    if isinstance(node.value, int):
        attr = ' type="int"'
        val = str(node.value)
    elif isinstance(node.value, float):
        attr = ' type="float"'
        val = str(node.value)
    else:
        attr = ""
        val = node.value
    file.write(f"{spaces}<{key}{attr}>{val}</{key}>\n")


def obj2xml(obj: Any, file_name: str) -> None:
    """Serialize a Python object to an XML file.

    Parameters
    ----------
    obj : Any
        Object whose ``__dict__`` will be serialized.
    file_name : str
        Output file path.
    """
    with open(file_name, "w") as f:
        f.write('<?xml version="1.0"?>\n<root>\n')
        write_node(f, obj, 2)
        f.write("</root>\n")


def dict2xml(xml_dict: dict, file_name: str) -> None:
    """Serialize a dictionary to an XML file.

    Parameters
    ----------
    xml_dict : dict
        Dictionary to serialize.
    file_name : str
        Output file path.
    """
    with open(file_name, "w") as f:
        f.write('<?xml version="1.0"?>\n<root>\n')
        write_node(f, xml_dict, 2)
        f.write("</root>\n")


def xml2obj(file_name: str) -> Any:
    """Deserialize an XML file (local or URL) to a Python object.

    Parameters
    ----------
    file_name : str
        Local file path or HTTP(S) URL.

    Returns
    -------
    Any
        Dynamically typed Python object tree.
    """
    if file_name.startswith("http"):
        with urllib.request.urlopen(file_name) as f:
            tree = ET.parse(f)
            xml_root = tree.getroot()
    else:
        xml_root = ET.parse(file_name).getroot()
    return _xml2py(xml_root)


def get_value(file_name: str, tag: str) -> Any:
    """Read a single value from an XML file by tag path.

    Parameters
    ----------
    file_name : str
        XML file path.
    tag : str
        XPath-style tag (e.g. ``"section/param"``).

    Returns
    -------
    Any
        The value, with type conversion applied if a ``type`` attribute exists.
    """
    xml_root = ET.parse(file_name).getroot()
    node = xml_root.findall(tag)[0]
    val = node.text
    if node.attrib and "type" in node.attrib:
        val = _convert_value(node.text, node.attrib["type"])
    return val


def _convert_value(text: str, type_name: str) -> Any:
    """Convert a string value based on a type name."""
    if type_name == "float":
        return float(text)
    elif type_name == "int":
        return int(text)
    elif type_name == "datetime":
        return datetime.strptime(text, "%Y%m%d %H%M%S")
    return text


def _xml2py(node: ET.Element) -> Any:
    """Recursively convert an XML element to a Python object."""
    name = node.tag
    pytype = type(name, (object,), {})
    pyobj = pytype()

    for attr in node.attrib:
        setattr(pyobj, attr, node.get(attr))

    if node.text and node.text.strip() not in ("", "\n"):
        setattr(pyobj, "text", node.text)
        setattr(pyobj, "value", node.text)
        if node.attrib and "type" in node.attrib:
            type_name = node.attrib["type"]
            if type_name == "float":
                lst = node.text.split(",")
                pyobj.value = (
                    float(node.text) if len(lst) == 1 else [float(s) for s in lst]
                )
            elif type_name == "int":
                if "," in node.text:
                    pyobj.value = [int(s) for s in node.text.split(",")]
                else:
                    pyobj.value = int(node.text)
            elif type_name == "datetime":
                pyobj.value = datetime.strptime(node.text, "%Y%m%d %H%M%S")

    for cn in node:
        if not hasattr(pyobj, cn.tag):
            setattr(pyobj, cn.tag, [])
        getattr(pyobj, cn.tag).append(_xml2py(cn))

    return pyobj

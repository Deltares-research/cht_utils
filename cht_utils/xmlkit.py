# -*- coding: utf-8 -*-
"""
Created on Wed May 12 08:59:03 2021

@author: ormondt
"""

import xml.etree.ElementTree as ET
from datetime import datetime
import urllib


class XMLObject:
    def __init__(self):
        # Empty object
        pass


def write_node(file, node, nspaces):

    if isinstance(node, dict):
        node_dict = node
        spaces = " " * nspaces
        for key in node_dict:
            for node2 in node_dict[key]:
                if not isinstance(node2, dict):
                    value = node
                    # End node
                    if type(node2.value) is int:
                        attrstr = ' type="int"'
                        valstr = str(node2.value)
                    elif type(node2.value) is float:
                        attrstr = ' type="float"'
                        valstr = str(node2.value)
                    else:
                        attrstr = ""
                        valstr = node2.value
                    str1 = "<" + key + attrstr + ">"
                    str2 = "</" + key + ">"
                    file.write(spaces + str1 + valstr + str2 + "\n")

                else:
                    file.write(spaces + "<" + key + ">\n")
                    write_node(file, node2, nspaces + 2)
                    file.write(spaces + "</" + key + ">\n")

    else:
        # Object
        node_dict = node.__dict__
        spaces = " " * nspaces
        for key in node_dict:
            for node2 in node_dict[key]:
                if hasattr(node2, "value"):
                    # End node
                    if type(node2.value) is int:
                        attrstr = ' type="int"'
                        valstr = str(node2.value)
                    elif type(node2.value) is float:
                        attrstr = ' type="float"'
                        valstr = str(node2.value)
                    else:
                        attrstr = ""
                        valstr = node2.value
                    str1 = "<" + key + attrstr + ">"
                    str2 = "</" + key + ">"
                    file.write(spaces + str1 + valstr + str2 + "\n")

                else:
                    file.write(spaces + "<" + key + ">\n")
                    write_node(file, node2, nspaces + 2)
                    file.write(spaces + "</" + key + ">\n")


def obj2xml(obj, file_name):

    file = open(file_name, "w")
    file.write('<?xml version="1.0"?>\n')
    file.write("<root>\n")
    write_node(file, obj, 2)
    file.write("</root>\n")
    file.close()


def dict2xml(xml_dict, file_name):

    file = open(file_name, "w")
    file.write('<?xml version="1.0"?>\n')
    file.write("<root>\n")
    write_node(file, xml_dict, 2)
    file.write("</root>\n")
    file.close()


def xml2py(node):

    name = node.tag

    pytype = type(name, (object,), {})
    pyobj = pytype()

    for attr in node.attrib.keys():
        setattr(pyobj, attr, node.get(attr))

    if node.text and node.text.strip() != "" and node.text.strip() != "\n":
        setattr(pyobj, "text", node.text)
        setattr(pyobj, "value", node.text)
        # Convert
        if node.attrib:
            if "type" in node.attrib.keys():
                if node.attrib["type"] == "float":
                    lst = node.text.split(",")
                    if len(lst) == 1:
                        pyobj.value = float(node.text)
                    else:
                        float_list = [float(s) for s in lst]
                        pyobj.value = float_list
                elif node.attrib["type"] == "int":
                    if "," in node.text:
                        pyobj.value = [int(s) for s in node.text.split(",")]
                    else:
                        pyobj.value = int(node.text)
                elif node.attrib["type"] == "datetime":
                    pyobj.value = datetime.strptime(node.text, "%Y%m%d %H%M%S")

    for cn in node:
        if not hasattr(pyobj, cn.tag):
            setattr(pyobj, cn.tag, [])
        getattr(pyobj, cn.tag).append(xml2py(cn))

    return pyobj


def get_value(file_name, tag):

    xml_root = ET.parse(file_name).getroot()
    node = xml_root.findall(tag)[0]
    val = node.text
    if node.attrib:
        if "type" in node.attrib.keys():
            if node.attrib["type"] == "float":
                val = float(node.text)
            elif node.attrib["type"] == "int":
                val = int(node.text)
            elif node.attrib["type"] == "datetime":
                val = datetime.strptime(node.text, "%Y%m%d %H%M%S")
    return val


def xml2obj(file_name):

    if file_name[0:4] == "http":
        with urllib.request.urlopen(file_name) as f:
            tree = ET.parse(f)
            xml_root = tree.getroot()
    else:
        xml_root = ET.parse(file_name).getroot()
    obj = xml2py(xml_root)

    return obj

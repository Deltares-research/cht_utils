"""Sphinx configuration for cht_utils."""

project = "cht_utils"
copyright = "2024, Deltares"
author = "Maarten van Ormondt"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
]

templates_path = ["_templates"]
exclude_patterns = ["_build"]

html_theme = "sphinx_rtd_theme"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
}

napoleon_numpy_docstring = True
autodoc_member_order = "bysource"

# Coastal Hazards Toolkit - Utils

Welcome to the GitHub page of the Deltares CHT Utils

For an editable package, install the package with e.g.:

cd d:\checkouts\github\cht_utils
pip install -e .

or something like:

pip install -e d:\checkouts\github\cht_utils






To build the cht_utils package, open Anaconda Powershell.

To make sure you have the latest version of build:

python -m pip install --upgrade build

and to build the package, e.g.: 

cd d:\checkouts\github\cht_utils

python -m build

Upload to Pypi with:

cd d:\checkouts\github\cht_utils

python -m pip install --upgrade twine

python -m twine upload dist/*

"""File and directory operations."""

# Core API
from cht_utils.fileops.fileops import copy as copy
from cht_utils.fileops.fileops import delete as delete
from cht_utils.fileops.fileops import exists as exists
from cht_utils.fileops.fileops import file_size as file_size
from cht_utils.fileops.fileops import find_replace as find_replace
from cht_utils.fileops.fileops import list_files as list_files
from cht_utils.fileops.fileops import list_folders as list_folders
from cht_utils.fileops.fileops import mkdir as mkdir
from cht_utils.fileops.fileops import move as move
from cht_utils.fileops.fileops import rename as rename
from cht_utils.fileops.fileops import touch as touch

# Backward-compatible aliases
from cht_utils.fileops.fileops import copy_file as copy_file
from cht_utils.fileops.fileops import delete_file as delete_file
from cht_utils.fileops.fileops import delete_folder as delete_folder
from cht_utils.fileops.fileops import findreplace as findreplace
from cht_utils.fileops.fileops import list_all_files as list_all_files
from cht_utils.fileops.fileops import list_files_recursive as list_files_recursive
from cht_utils.fileops.fileops import move_file as move_file
from cht_utils.fileops.fileops import rm as rm
from cht_utils.fileops.fileops import rmdir as rmdir

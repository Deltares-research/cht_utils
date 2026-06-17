"""File and directory operations for copying, moving, deleting, and listing paths."""

from __future__ import annotations

import glob
import logging
import os
import shutil

logger = logging.getLogger(__name__)


def move_file(src: str, dst: str) -> None:
    """Move files matching a glob pattern to a destination directory.

    Parameters
    ----------
    src : str
        Glob pattern for source files.
    dst : str
        Destination directory path.
    """
    for full_file_name in glob.glob(src):
        src_name = os.path.basename(full_file_name)
        if os.path.exists(os.path.join(dst, src_name)):
            try:
                os.remove(os.path.join(dst, src_name))
            except Exception:
                logger.error(f"Could not remove file {os.path.join(dst, src_name)}")
        try:
            shutil.move(full_file_name, dst)
        except Exception:
            logger.error(f"Could not move file {full_file_name}")


def copy_file(src: str, dst: str) -> None:
    """Copy files matching a glob pattern to a destination directory.

    Parameters
    ----------
    src : str
        Glob pattern for source files.
    dst : str
        Destination directory path.
    """
    for full_file_name in glob.glob(src):
        src_name = os.path.basename(full_file_name)
        if os.path.exists(os.path.join(dst, src_name)):
            os.remove(os.path.join(dst, src_name))
        if os.path.isdir(full_file_name):
            dstf = os.path.join(dst, os.path.basename(full_file_name))
            shutil.copytree(full_file_name, dstf)
        else:
            shutil.copy(full_file_name, dst)


def delete_file(src: str) -> None:
    """Delete files matching a glob pattern.

    Parameters
    ----------
    src : str
        Glob pattern for files to delete.
    """
    for file_name in glob.glob(src):
        try:
            os.remove(src)
        except Exception:
            logger.error(f"Could not delete {src}")


def rm(src: str) -> None:
    """Remove a single file.

    Parameters
    ----------
    src : str
        Path to the file to remove.
    """
    os.remove(src)


def mkdir(path: str) -> None:
    """Create a directory (and parents) if it does not already exist.

    Parameters
    ----------
    path : str
        Directory path to create.
    """
    if not os.path.exists(path):
        os.makedirs(path)


def list_files(src: str, full_path: bool = True) -> list[str]:
    """List files matching a glob pattern or directory listing.

    Parameters
    ----------
    src : str
        Glob pattern (when *full_path* is True) or directory path (when False).
    full_path : bool
        If True, use glob and return full paths. If False, use ``os.listdir``.

    Returns
    -------
    list[str]
        Sorted list of file paths or names.
    """
    if full_path:
        file_list = []
        full_list = glob.glob(src)
        for item in full_list:
            if os.path.isfile(item):
                file_list.append(item)
    else:
        file_list = os.listdir(src)

    return sorted(file_list)


def list_folders(src: str, basename: bool = False) -> list[str]:
    """List directories matching a glob pattern.

    Parameters
    ----------
    src : str
        Glob pattern.
    basename : bool
        If True, return only base names instead of full paths.

    Returns
    -------
    list[str]
        Sorted list of directory paths (or base names).
    """
    folder_list = []
    full_list = glob.glob(src)
    for item in full_list:
        if os.path.isdir(item):
            if basename:
                folder_list.append(os.path.basename(item))
            else:
                folder_list.append(item)

    return sorted(folder_list)


def delete_folder(src: str) -> None:
    """Recursively delete a directory tree.

    Parameters
    ----------
    src : str
        Path to the directory to delete.
    """
    try:
        shutil.rmtree(src, ignore_errors=False, onerror=None)
    except Exception:
        logger.error(f"Could not delete folder {src}")


def rmdir(src: str) -> None:
    """Recursively delete a directory tree if it exists.

    Parameters
    ----------
    src : str
        Path to the directory to delete.
    """
    try:
        if os.path.exists(src):
            shutil.rmtree(src, ignore_errors=False, onerror=None)
    except Exception:
        logger.error(f"Could not delete folder {src}")


def exists(src: str) -> bool:
    """Check whether a path exists.

    Parameters
    ----------
    src : str
        Path to check.

    Returns
    -------
    bool
        True if the path exists.
    """
    return os.path.exists(src)

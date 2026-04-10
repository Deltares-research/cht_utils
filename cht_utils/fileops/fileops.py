"""Common file and directory operations using pathlib and shutil."""

import shutil
from pathlib import Path
from typing import List, Union

PathLike = Union[str, Path]


def move(src: PathLike, dst: PathLike) -> None:
    """Move files to a destination.

    If *src* contains glob wildcards, all matching files are moved into *dst*
    (which must be a directory). Without wildcards, *dst* may be a directory
    or a new file name (rename).

    Parameters
    ----------
    src : str or Path
        Source glob pattern (e.g. ``"data/*.nc"``) or single file/directory.
    dst : str or Path
        Destination directory or new name.
    """
    src = Path(src)
    dst = Path(dst)
    if _has_glob(src):
        for f in src.parent.glob(src.name):
            target = dst / f.name
            if target.exists():
                if target.is_dir():
                    shutil.rmtree(target)
                else:
                    target.unlink()
            shutil.move(str(f), str(dst))
    else:
        if src.exists():
            # Check if dst is an existing directory; if so, move src into it
            # if dst.is_dir():
            #     dst = dst / src.name
            shutil.move(str(src), str(dst))


def copy(src: PathLike, dst: PathLike) -> None:
    """Copy files or directories to a destination.

    If *src* contains glob wildcards, all matching files are copied into *dst*
    (which must be a directory). Without wildcards, *dst* may be a directory
    or a new file name.

    Parameters
    ----------
    src : str or Path
        Source glob pattern or single file/directory.
    dst : str or Path
        Destination directory or new name.
    """
    src = Path(src)
    dst = Path(dst)
    if _has_glob(src):
        for f in src.parent.glob(src.name):
            _copy_single(f, dst / f.name)
    else:
        if not src.exists():
            print(f"{src} does not exist")
            return
        if dst.is_dir():
            _copy_single(src, dst / src.name)
        else:
            _copy_single(src, dst)


def delete(src: Union[PathLike, List[PathLike]]) -> None:
    """Delete files or directories matching glob patterns.

    Handles both files (``unlink``) and directories (``rmtree``).
    Silently skips paths that don't exist.

    Parameters
    ----------
    src : str, Path, or list thereof
        Glob pattern(s) or explicit path(s) to delete.
    """
    patterns = src if isinstance(src, list) else [src]
    for pattern in patterns:
        p = Path(pattern)
        if _has_glob(p):
            matches = list(p.parent.glob(p.name))
        else:
            matches = [p] if p.exists() else []
        for f in matches:
            try:
                if f.is_dir():
                    shutil.rmtree(f)
                else:
                    f.unlink()
            except OSError:
                print(f"Could not delete {f}")


def rename(src: PathLike, dst: PathLike) -> Path:
    """Rename a file or directory.

    Works across drives on Windows (falls back to copy + delete).

    Parameters
    ----------
    src : str or Path
        Source path.
    dst : str or Path
        New name or path.

    Returns
    -------
    Path
        The new path.
    """
    src = Path(src)
    dst = Path(dst)
    try:
        return src.rename(dst)
    except OSError:
        # Cross-drive rename on Windows
        if src.is_dir():
            shutil.copytree(src, dst)
            shutil.rmtree(src)
        else:
            shutil.copy2(src, dst)
            src.unlink()
        return dst


def mkdir(path: PathLike) -> Path:
    """Create a directory (including parents) if it does not exist.

    Parameters
    ----------
    path : str or Path
        Directory path to create.

    Returns
    -------
    Path
        The created directory path.
    """
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def touch(path: PathLike) -> Path:
    """Create an empty file or update its modification timestamp.

    Parent directories are created if needed.

    Parameters
    ----------
    path : str or Path
        File path.

    Returns
    -------
    Path
        The touched file path.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.touch()
    return p


def file_size(path: PathLike) -> int:
    """Return the size of a file in bytes.

    Parameters
    ----------
    path : str or Path
        File path.

    Returns
    -------
    int
        File size in bytes.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    """
    return Path(path).stat().st_size


def list_files(
    src: PathLike,
    pattern: str = "*",
    recursive: bool = False,
    full_path: bool = True,
) -> List[str]:
    """List files in a directory, optionally filtered by glob pattern.

    Parameters
    ----------
    src : str or Path
        Directory path, or glob pattern for backward compatibility
        (e.g. ``"data/*.nc"``).
    pattern : str
        Glob pattern to filter files. Defaults to ``"*"`` (all files).
        Ignored when *src* itself contains wildcards.
    recursive : bool
        If ``True``, search subdirectories recursively.
    full_path : bool
        If ``True``, return full paths. If ``False``, return basenames.

    Returns
    -------
    List[str]
        Sorted list of file paths or names.
    """

    # Collapse consecutive wildcards — pathlib rejects "**" mixed with text
    while "**" in pattern:
        pattern = pattern.replace("**", "*")
    p = Path(src)
    if _has_glob(p):
        # Legacy: src is a glob pattern like "data/*.nc"
        name = p.name
        while "**" in name:
            name = name.replace("**", "*")
        files = list(p.parent.glob(name))
    elif recursive:
        files = list(p.rglob(pattern))
    else:
        files = list(p.glob(pattern))
    files = [f for f in files if f.is_file()]
    if full_path:
        return sorted(str(f) for f in files)
    return sorted(f.name for f in files)


def list_folders(
    src: PathLike, pattern: str = "*", basename: bool = False
) -> List[str]:
    """List subdirectories, optionally filtered by glob pattern.

    Parameters
    ----------
    src : str or Path
        Directory path, or glob pattern for backward compatibility.
    pattern : str
        Glob pattern to filter folders. Defaults to ``"*"``.
        Ignored when *src* itself contains wildcards.
    basename : bool
        If ``True``, return only folder names instead of full paths.

    Returns
    -------
    List[str]
        Sorted list of folder paths or names.
    """
    p = Path(src)
    if _has_glob(p):
        folders = list(p.parent.glob(p.name))
    else:
        folders = list(p.glob(pattern))
    folders = [f for f in folders if f.is_dir()]
    if basename:
        return sorted(f.name for f in folders)
    return sorted(str(f) for f in folders)


def exists(src: PathLike) -> bool:
    """Check if a file or directory exists.

    Parameters
    ----------
    src : str or Path
        Path to check.

    Returns
    -------
    bool
    """
    return Path(src).exists()


def find_replace(file_name: PathLike, old: str, new: str) -> None:
    """Replace all occurrences of a string in a file.

    Parameters
    ----------
    file_name : str or Path
        File to modify in-place.
    old : str
        String to search for.
    new : str
        Replacement string.
    """
    p = Path(file_name)
    p.write_text(p.read_text().replace(old, new))


# ── Internal helpers ──────────────────────────────────────────────────


def _has_glob(p: Path) -> bool:
    """Check if a path contains glob wildcard characters."""
    return any(c in str(p) for c in ("*", "?", "["))


def _copy_single(src: Path, dst: Path) -> None:
    """Copy a single file or directory to a destination path."""
    if dst.exists():
        if dst.is_dir():
            shutil.rmtree(dst)
        else:
            dst.unlink()
    if src.is_dir():
        shutil.copytree(src, dst)
    else:
        shutil.copy2(src, dst)


# ── Backward-compatible aliases ───────────────────────────────────────

move_file = move
copy_file = copy
delete_file = delete
delete_folder = delete
rm = delete
rmdir = delete
findreplace = find_replace


def list_all_files(src: PathLike) -> List[str]:
    """Backward-compatible alias for ``list_files(src, recursive=True)``."""
    return list_files(src, recursive=True)


def list_files_recursive(src: PathLike, pattern: str = "*") -> List[str]:
    """Backward-compatible alias for ``list_files(src, pattern, recursive=True)``."""
    return list_files(src, pattern=pattern, recursive=True)

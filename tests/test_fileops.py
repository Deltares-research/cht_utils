"""Tests for cht_utils.fileops."""

import pytest

from cht_utils import fileops as fo


class TestMkdir:
    def test_creates_directory(self, tmp_path):
        d = tmp_path / "a" / "b" / "c"
        fo.mkdir(str(d))
        assert d.exists()

    def test_existing_directory(self, tmp_path):
        fo.mkdir(str(tmp_path))  # should not raise

    def test_returns_path(self, tmp_path):
        result = fo.mkdir(tmp_path / "new")
        assert result.exists()


class TestTouch:
    def test_creates_empty_file(self, tmp_path):
        f = tmp_path / "new.txt"
        fo.touch(str(f))
        assert f.exists()
        assert f.stat().st_size == 0

    def test_creates_parent_dirs(self, tmp_path):
        f = tmp_path / "a" / "b" / "file.txt"
        fo.touch(str(f))
        assert f.exists()

    def test_returns_path(self, tmp_path):
        result = fo.touch(tmp_path / "x.txt")
        assert result.exists()


class TestFileSize:
    def test_empty_file(self, tmp_path):
        f = tmp_path / "empty.txt"
        f.touch()
        assert fo.file_size(str(f)) == 0

    def test_file_with_content(self, tmp_path):
        f = tmp_path / "data.txt"
        f.write_text("hello")
        assert fo.file_size(str(f)) == 5

    def test_nonexistent_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            fo.file_size(str(tmp_path / "nope.txt"))


class TestCopy:
    def test_copy_single_file(self, tmp_path):
        src = tmp_path / "src.txt"
        dst = tmp_path / "dst.txt"
        src.write_text("hello")
        fo.copy(str(src), str(dst))
        assert dst.read_text() == "hello"

    def test_copy_to_directory(self, tmp_path):
        src = tmp_path / "src.txt"
        dst_dir = tmp_path / "out"
        dst_dir.mkdir()
        src.write_text("data")
        fo.copy(str(src), str(dst_dir))
        assert (dst_dir / "src.txt").read_text() == "data"

    def test_copy_wildcard(self, tmp_path):
        (tmp_path / "a.nc").write_text("a")
        (tmp_path / "b.nc").write_text("b")
        (tmp_path / "c.txt").write_text("c")
        dst = tmp_path / "out"
        dst.mkdir()
        fo.copy(str(tmp_path / "*.nc"), str(dst))
        assert (dst / "a.nc").exists()
        assert (dst / "b.nc").exists()
        assert not (dst / "c.txt").exists()

    def test_copy_directory(self, tmp_path):
        src = tmp_path / "mydir"
        src.mkdir()
        (src / "file.txt").write_text("inside")
        dst = tmp_path / "out"
        dst.mkdir()
        fo.copy(str(src), str(dst))
        assert (dst / "mydir" / "file.txt").read_text() == "inside"


class TestMove:
    def test_move_file(self, tmp_path):
        src = tmp_path / "src.txt"
        dst_dir = tmp_path / "out"
        dst_dir.mkdir()
        src.write_text("data")
        fo.move(str(src), str(dst_dir))
        assert not src.exists()
        assert (dst_dir / "src.txt").read_text() == "data"

    def test_move_rename(self, tmp_path):
        src = tmp_path / "old.txt"
        dst = tmp_path / "new.txt"
        src.write_text("data")
        fo.move(str(src), str(dst))
        assert not src.exists()
        assert dst.read_text() == "data"

    def test_move_wildcard(self, tmp_path):
        (tmp_path / "a.nc").write_text("a")
        (tmp_path / "b.nc").write_text("b")
        dst = tmp_path / "out"
        dst.mkdir()
        fo.move(str(tmp_path / "*.nc"), str(dst))
        assert not (tmp_path / "a.nc").exists()
        assert (dst / "a.nc").exists()
        assert (dst / "b.nc").exists()


class TestDelete:
    def test_delete_file(self, tmp_path):
        f = tmp_path / "file.txt"
        f.write_text("data")
        fo.delete(str(f))
        assert not f.exists()

    def test_delete_directory(self, tmp_path):
        d = tmp_path / "mydir"
        d.mkdir()
        (d / "file.txt").write_text("data")
        fo.delete(str(d))
        assert not d.exists()

    def test_delete_wildcard(self, tmp_path):
        (tmp_path / "a.tmp").write_text("a")
        (tmp_path / "b.tmp").write_text("b")
        (tmp_path / "keep.txt").write_text("keep")
        fo.delete(str(tmp_path / "*.tmp"))
        assert not (tmp_path / "a.tmp").exists()
        assert not (tmp_path / "b.tmp").exists()
        assert (tmp_path / "keep.txt").exists()

    def test_delete_list(self, tmp_path):
        f1 = tmp_path / "a.txt"
        f2 = tmp_path / "b.txt"
        f1.write_text("a")
        f2.write_text("b")
        fo.delete([str(f1), str(f2)])
        assert not f1.exists()
        assert not f2.exists()

    def test_delete_nonexistent_silent(self, tmp_path):
        fo.delete(str(tmp_path / "nope.txt"))  # should not raise


class TestRename:
    def test_rename_file(self, tmp_path):
        src = tmp_path / "old.txt"
        dst = tmp_path / "new.txt"
        src.write_text("data")
        result = fo.rename(str(src), str(dst))
        assert not src.exists()
        assert dst.read_text() == "data"


class TestListFiles:
    def test_list_all_files(self, tmp_path):
        (tmp_path / "a.txt").touch()
        (tmp_path / "b.nc").touch()
        (tmp_path / "subdir").mkdir()
        result = fo.list_files(str(tmp_path))
        assert len(result) == 2  # excludes subdir

    def test_list_with_pattern(self, tmp_path):
        (tmp_path / "a.txt").touch()
        (tmp_path / "b.nc").touch()
        result = fo.list_files(str(tmp_path), "*.nc")
        assert len(result) == 1
        assert "b.nc" in result[0]

    def test_list_recursive(self, tmp_path):
        (tmp_path / "a.nc").touch()
        sub = tmp_path / "sub"
        sub.mkdir()
        (sub / "b.nc").touch()
        result = fo.list_files(str(tmp_path), "*.nc", recursive=True)
        assert len(result) == 2

    def test_list_basenames(self, tmp_path):
        (tmp_path / "a.txt").touch()
        (tmp_path / "b.txt").touch()
        result = fo.list_files(str(tmp_path), full_path=False)
        assert result == ["a.txt", "b.txt"]

    def test_legacy_glob_pattern(self, tmp_path):
        (tmp_path / "a.nc").touch()
        (tmp_path / "b.txt").touch()
        result = fo.list_files(str(tmp_path / "*.nc"))
        assert len(result) == 1

    def test_sorted_output(self, tmp_path):
        (tmp_path / "c.txt").touch()
        (tmp_path / "a.txt").touch()
        (tmp_path / "b.txt").touch()
        result = fo.list_files(str(tmp_path), full_path=False)
        assert result == ["a.txt", "b.txt", "c.txt"]


class TestListFolders:
    def test_list_folders(self, tmp_path):
        (tmp_path / "dir1").mkdir()
        (tmp_path / "dir2").mkdir()
        (tmp_path / "file.txt").touch()
        result = fo.list_folders(str(tmp_path))
        assert len(result) == 2

    def test_list_folders_basename(self, tmp_path):
        (tmp_path / "alpha").mkdir()
        (tmp_path / "beta").mkdir()
        result = fo.list_folders(str(tmp_path), basename=True)
        assert result == ["alpha", "beta"]

    def test_list_folders_with_pattern(self, tmp_path):
        (tmp_path / "data_01").mkdir()
        (tmp_path / "data_02").mkdir()
        (tmp_path / "other").mkdir()
        result = fo.list_folders(str(tmp_path), pattern="data_*", basename=True)
        assert result == ["data_01", "data_02"]

    def test_legacy_glob_pattern(self, tmp_path):
        (tmp_path / "dir1").mkdir()
        (tmp_path / "dir2").mkdir()
        result = fo.list_folders(str(tmp_path / "*"), basename=True)
        assert len(result) == 2


class TestExists:
    def test_file_exists(self, tmp_path):
        f = tmp_path / "file.txt"
        f.touch()
        assert fo.exists(str(f))

    def test_dir_exists(self, tmp_path):
        assert fo.exists(str(tmp_path))

    def test_nonexistent(self, tmp_path):
        assert not fo.exists(str(tmp_path / "nope"))


class TestFindReplace:
    def test_replaces_string(self, tmp_path):
        f = tmp_path / "file.txt"
        f.write_text("hello world, hello earth")
        fo.find_replace(str(f), "hello", "bye")
        assert f.read_text() == "bye world, bye earth"

    def test_no_match(self, tmp_path):
        f = tmp_path / "file.txt"
        f.write_text("nothing here")
        fo.find_replace(str(f), "xyz", "abc")
        assert f.read_text() == "nothing here"


class TestBackwardCompat:
    """Verify old function names still work."""

    def test_move_file(self, tmp_path):
        src = tmp_path / "src.txt"
        dst = tmp_path / "out"
        dst.mkdir()
        src.write_text("data")
        fo.move_file(str(src), str(dst))
        assert (dst / "src.txt").exists()

    def test_copy_file(self, tmp_path):
        src = tmp_path / "src.txt"
        dst = tmp_path / "out"
        dst.mkdir()
        src.write_text("data")
        fo.copy_file(str(src), str(dst))
        assert (dst / "src.txt").exists()

    def test_delete_file(self, tmp_path):
        f = tmp_path / "file.txt"
        f.touch()
        fo.delete_file(str(f))
        assert not f.exists()

    def test_delete_folder(self, tmp_path):
        d = tmp_path / "dir"
        d.mkdir()
        fo.delete_folder(str(d))
        assert not d.exists()

    def test_rm(self, tmp_path):
        f = tmp_path / "file.txt"
        f.touch()
        fo.rm(str(f))
        assert not f.exists()

    def test_rmdir(self, tmp_path):
        d = tmp_path / "dir"
        d.mkdir()
        fo.rmdir(str(d))
        assert not d.exists()

    def test_list_all_files(self, tmp_path):
        (tmp_path / "a.txt").touch()
        sub = tmp_path / "sub"
        sub.mkdir()
        (sub / "b.txt").touch()
        result = fo.list_all_files(str(tmp_path))
        assert len(result) == 2

    def test_findreplace(self, tmp_path):
        f = tmp_path / "file.txt"
        f.write_text("old text")
        fo.findreplace(str(f), "old", "new")
        assert f.read_text() == "new text"

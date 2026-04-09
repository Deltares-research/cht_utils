"""Tests for cht_utils.fileio subpackage."""

from cht_utils.fileio.deltares_ini import IniStruct
from cht_utils.fileio.json_js import read_json_js, write_csv_js, write_json_js
from cht_utils.fileio.yaml import dict2yaml, yaml2dict


class TestYaml:
    def test_roundtrip(self, tmp_path):
        f = tmp_path / "test.yml"
        data = {"name": "test", "values": [1, 2, 3], "nested": {"a": 1.5}}
        dict2yaml(str(f), data)
        result = yaml2dict(str(f))
        assert result["name"] == "test"
        assert result["values"] == [1, 2, 3]
        assert result["nested"]["a"] == 1.5

    def test_sort_keys(self, tmp_path):
        f = tmp_path / "sorted.yml"
        data = {"z": 1, "a": 2}
        dict2yaml(str(f), data, sort_keys=True)
        content = f.read_text()
        assert content.index("a:") < content.index("z:")


class TestJsonJs:
    def test_dict_roundtrip(self, tmp_path):
        f = tmp_path / "data.js"
        data = {"key": "value", "num": 42}
        write_json_js(str(f), data, "var data =")
        result = read_json_js(str(f))
        assert result["key"] == "value"
        assert result["num"] == 42

    def test_list_roundtrip(self, tmp_path):
        f = tmp_path / "list.js"
        data = [{"a": 1}, {"b": 2}]
        write_json_js(str(f), data, "var items =")
        result = read_json_js(str(f))
        assert len(result) == 2
        assert result[0]["a"] == 1

    def test_csv_js(self, tmp_path):
        f = tmp_path / "sub" / "data.js"
        write_csv_js(str(f), "a,b\n1,2\n", "var csv = `")
        assert f.exists()
        content = f.read_text()
        assert content.startswith("var csv = `")
        assert "a,b" in content


class TestDeltaresIni:
    def _write_ini(self, tmp_path):
        f = tmp_path / "test.ini"
        f.write_text(
            "[General]\n"
            "   name                 = MyModel\n"
            "   version              = 1.0          # model version\n"
            "\n"
            "[Output]\n"
            "   format               = netcdf\n"
            "   interval             = 3600\n"
        )
        return str(f)

    def test_read_sections(self, tmp_path):
        ini = IniStruct(self._write_ini(tmp_path))
        assert len(ini.section) == 2
        assert ini.section[0].name == "General"
        assert ini.section[1].name == "Output"

    def test_read_keywords(self, tmp_path):
        ini = IniStruct(self._write_ini(tmp_path))
        assert ini.get_value("General", "name") == "MyModel"
        assert ini.get_value("General", "version") == "1.0"
        assert ini.get_value("Output", "format") == "netcdf"

    def test_read_comment(self, tmp_path):
        ini = IniStruct(self._write_ini(tmp_path))
        kw = ini.section[0].keyword[1]
        assert kw.comment == "model version"

    def test_set_existing_value(self, tmp_path):
        ini = IniStruct(self._write_ini(tmp_path))
        ini.set_value("General", "name", "NewName")
        assert ini.get_value("General", "name") == "NewName"

    def test_set_new_keyword(self, tmp_path):
        ini = IniStruct(self._write_ini(tmp_path))
        ini.set_value("General", "author", "test", "added by test")
        assert ini.get_value("General", "author") == "test"

    def test_set_new_section(self, tmp_path):
        ini = IniStruct(self._write_ini(tmp_path))
        ini.set_value("NewSection", "key", "val")
        assert ini.get_value("NewSection", "key") == "val"

    def test_write_roundtrip(self, tmp_path):
        ini = IniStruct(self._write_ini(tmp_path))
        out = tmp_path / "out.ini"
        ini.write(str(out))
        ini2 = IniStruct(str(out))
        assert ini2.get_value("General", "name") == "MyModel"
        assert ini2.get_value("Output", "interval") == "3600"

    def test_get_section(self, tmp_path):
        ini = IniStruct(self._write_ini(tmp_path))
        sec = ini.get_section("Output")
        assert sec.name == "Output"

    def test_get_value_missing(self, tmp_path):
        ini = IniStruct(self._write_ini(tmp_path))
        assert ini.get_value("General", "nonexistent") is None
        assert ini.get_value("NonExistent", "key") is None

    def test_section_get_value(self, tmp_path):
        ini = IniStruct(self._write_ini(tmp_path))
        val = ini.section[0].get_value("name")
        assert val == "MyModel"


class TestImports:
    """Verify subpackage re-exports work."""

    def test_fileio_imports(self):
        pass

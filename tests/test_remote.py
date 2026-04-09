"""Tests for cht_utils.remote subpackage.

These tests only check import and instantiation logic, not actual
network operations (which require credentials and connectivity).
"""

import pytest


class TestS3SessionImport:
    def test_import(self):
        from cht_utils.remote.s3 import S3Session
        assert S3Session is not None

    def test_import_from_package(self):
        from cht_utils.remote import S3Session
        assert S3Session is not None


class TestSSHSessionImport:
    def test_import(self):
        from cht_utils.remote.sftp import SSHSession
        assert SSHSession is not None

    def test_import_from_package(self):
        from cht_utils.remote import SSHSession
        assert SSHSession is not None

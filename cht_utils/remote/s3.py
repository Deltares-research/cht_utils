"""AWS S3 client wrapper for file and folder operations."""

import os
import pathlib
import tarfile
from multiprocessing.pool import ThreadPool
from typing import List, Optional

import boto3
from botocore.exceptions import ClientError


class S3Session:
    """AWS S3 session with upload, download, and management operations.

    Parameters
    ----------
    access_key : str
        AWS access key ID.
    secret_key : str
        AWS secret access key.
    region : str
        AWS region name.
    """

    def __init__(self, access_key: str, secret_key: str, region: str) -> None:
        self.ready = True
        session = boto3.Session(
            aws_access_key_id=access_key,
            aws_secret_access_key=secret_key,
            region_name=region,
        )
        self.s3_client = session.client("s3")

    def upload_file(
        self, bucket_name: str, file: str, s3_folder: str, quiet: bool = True
    ) -> None:
        """Upload a local file to S3.

        Parameters
        ----------
        bucket_name : str
            S3 bucket name.
        file : str
            Local file path.
        s3_folder : str
            Remote folder prefix.
        quiet : bool
            Suppress progress output.
        """
        s3_key = os.path.join(s3_folder, os.path.basename(file)).replace("\\", "/")
        self.s3_client.upload_file(file, bucket_name, s3_key)
        if not quiet:
            print(f"Uploaded {os.path.basename(file)}")

    def download_file(
        self,
        bucket_name: str,
        s3_folder: str,
        file: str,
        local_folder: str,
        quiet: bool = True,
    ) -> None:
        """Download a file from S3.

        Parameters
        ----------
        bucket_name : str
            S3 bucket name.
        s3_folder : str
            Remote folder prefix.
        file : str
            File name to download.
        local_folder : str
            Local destination directory.
        quiet : bool
            Suppress progress output.
        """
        s3_key = os.path.join(s3_folder, os.path.basename(file)).replace("\\", "/")
        local_path = os.path.join(local_folder, os.path.basename(file))
        self.s3_client.download_file(bucket_name, s3_key, local_path)
        if not quiet:
            print(f"Downloaded {os.path.basename(file)}")

    def download_files(
        self,
        bucket_name: str,
        key_list: List[str],
        file_list: List[str],
        quiet: bool = True,
    ) -> None:
        """Download multiple files in parallel.

        Parameters
        ----------
        bucket_name : str
            S3 bucket name.
        key_list : List[str]
            List of S3 keys.
        file_list : List[str]
            Corresponding local file paths.
        quiet : bool
            Suppress progress output.
        """
        pool = ThreadPool()
        pool.starmap(
            self.download_file_parallel,
            [(bucket_name, key, file, quiet) for key, file in zip(key_list, file_list)],
        )

    def download_file_parallel(
        self, bucket_name: str, s3_key: str, file: str, quiet: bool = True
    ) -> None:
        """Download a single file (used for parallel execution).

        Parameters
        ----------
        bucket_name : str
            S3 bucket name.
        s3_key : str
            Full S3 key.
        file : str
            Local file path.
        quiet : bool
            Suppress progress output.
        """
        self.s3_client.download_file(bucket_name, s3_key, file)
        if not quiet:
            print(f"Downloaded {os.path.basename(file)}")

    def delete_file(
        self, bucket_name: str, s3_folder: str, file: str, quiet: bool = True
    ) -> None:
        """Delete a file from S3.

        Parameters
        ----------
        bucket_name : str
            S3 bucket name.
        s3_folder : str
            Remote folder prefix.
        file : str
            File name to delete.
        quiet : bool
            Suppress progress output.
        """
        s3_key = os.path.join(s3_folder, os.path.basename(file)).replace("\\", "/")
        self.s3_client.delete_object(Bucket=bucket_name, Key=s3_key)
        if not quiet:
            print(f"Deleted {os.path.basename(file)}")

    def make_folder(
        self, bucket_name: str, s3_folder: str, quiet: bool = True
    ) -> None:
        """Create a folder (prefix) in S3.

        Parameters
        ----------
        bucket_name : str
            S3 bucket name.
        s3_folder : str
            Folder prefix to create.
        quiet : bool
            Suppress progress output.
        """
        self.s3_client.put_object(Bucket=bucket_name, Key=f"{s3_folder}/")
        if not quiet:
            print(f"Made folder: {s3_folder}")

    def upload_folder(
        self,
        bucket_name: str,
        local_folder: str,
        s3_folder: str,
        parallel: bool = True,
        quiet: bool = True,
    ) -> None:
        """Recursively upload a local directory to S3.

        Parameters
        ----------
        bucket_name : str
            S3 bucket name.
        local_folder : str
            Local directory to upload.
        s3_folder : str
            Remote folder prefix.
        parallel : bool
            Use parallel uploads.
        quiet : bool
            Suppress progress output.
        """
        local_folder = local_folder.replace("\\", "/")
        flist = _list_all_files(local_folder)
        args = [
            (f, local_folder, s3_folder, bucket_name, self.s3_client, quiet)
            for f in flist
        ]
        if parallel:
            pool = ThreadPool()
            pool.starmap(_upload_file, args)
        else:
            for a in args:
                _upload_file(*a)

    def download_folder(
        self,
        bucket_name: str,
        s3_folder: str,
        local_folder: str,
        parallel: bool = True,
        quiet: bool = True,
    ) -> None:
        """Download all files under an S3 prefix to a local directory.

        Parameters
        ----------
        bucket_name : str
            S3 bucket name.
        s3_folder : str
            Remote folder prefix.
        local_folder : str
            Local destination directory.
        parallel : bool
            Currently unused (sequential download).
        quiet : bool
            Suppress progress output.
        """
        os.makedirs(local_folder, exist_ok=True)
        objects = self.s3_client.list_objects(Bucket=bucket_name, Prefix=s3_folder)
        for obj in objects.get("Contents", []):
            s3_key = obj["Key"]
            local_path = os.path.join(local_folder, os.path.basename(s3_key))
            self.s3_client.download_file(bucket_name, s3_key, local_path)
            if not quiet:
                print(f"Downloaded {os.path.basename(s3_key)}")

    def delete_folder(self, bucket_name: str, s3_folder: str) -> None:
        """Delete all objects under an S3 prefix.

        Parameters
        ----------
        bucket_name : str
            S3 bucket name.
        s3_folder : str
            Remote folder prefix.
        """
        if not s3_folder.endswith("/"):
            s3_folder += "/"
        objects = self.s3_client.list_objects(Bucket=bucket_name, Prefix=s3_folder)
        for obj in objects.get("Contents", []):
            self.s3_client.delete_object(Bucket=bucket_name, Key=obj["Key"])

    def list_folders(self, bucket_name: str, s3_folder: str) -> List[str]:
        """List immediate subfolders under an S3 prefix.

        Parameters
        ----------
        bucket_name : str
            S3 bucket name.
        s3_folder : str
            Remote folder prefix.

        Returns
        -------
        List[str]
            Subfolder names.
        """
        if not s3_folder.endswith("/"):
            s3_folder += "/"
        folders = []
        paginator = self.s3_client.get_paginator("list_objects_v2")
        for page in paginator.paginate(
            Bucket=bucket_name, Prefix=s3_folder, Delimiter="/"
        ):
            for subfolder in page.get("CommonPrefixes", []):
                folders.append(subfolder["Prefix"].rstrip("/").split("/")[-1])
        return folders

    def list_files(self, bucket_name: str, s3_folder: str) -> List[str]:
        """List all file keys under an S3 prefix.

        Parameters
        ----------
        bucket_name : str
            S3 bucket name.
        s3_folder : str
            Remote folder prefix.

        Returns
        -------
        List[str]
            Full S3 keys.
        """
        all_files = []
        paginator = self.s3_client.get_paginator("list_objects_v2")
        for page in paginator.paginate(Bucket=bucket_name, Prefix=s3_folder):
            for obj in page.get("Contents", []):
                all_files.append(obj["Key"])
        return all_files

    def download_and_extract_tgz(
        self, bucket_name: str, s3_folder: str, local_folder: str
    ) -> None:
        """Download and extract a ``.tgz`` archive from S3.

        Parameters
        ----------
        bucket_name : str
            S3 bucket name.
        s3_folder : str
            S3 key of the ``.tgz`` file.
        local_folder : str
            Local directory to extract into.
        """
        local_tgz_path = os.path.join("/tmp", os.path.basename(s3_folder))
        self.s3_client.download_file(bucket_name, s3_folder, local_tgz_path)
        with tarfile.open(local_tgz_path, "r:gz") as tar:
            tar.extractall(path=local_folder)
        os.remove(local_tgz_path)

    def check_file_exists(self, bucket_name: str, s3_key: str) -> bool:
        """Check if a file exists in S3.

        Parameters
        ----------
        bucket_name : str
            S3 bucket name.
        s3_key : str
            Full S3 key.

        Returns
        -------
        bool
        """
        try:
            self.s3_client.head_object(Bucket=bucket_name, Key=s3_key)
            return True
        except ClientError as e:
            if e.response["Error"]["Code"] == "404":
                return False
            raise

    def check_folder_exists(self, bucket_name: str, s3_key: str) -> bool:
        """Check if a folder (prefix) exists in S3.

        Parameters
        ----------
        bucket_name : str
            S3 bucket name.
        s3_key : str
            S3 prefix to check.

        Returns
        -------
        bool
        """
        response = self.s3_client.list_objects_v2(
            Bucket=bucket_name, Prefix=s3_key, Delimiter="/"
        )
        return "CommonPrefixes" in response


def _list_all_files(src: str) -> List[str]:
    """Recursively list all files in a directory."""
    return [str(f) for f in pathlib.Path(src).rglob("*") if f.is_file()]


def _upload_file(
    file: str,
    local_folder: str,
    s3_folder: str,
    bucket_name: str,
    s3_client: object,
    quiet: bool,
) -> None:
    """Upload a single file, preserving relative path structure."""
    relative = file.replace("\\", "/").replace(local_folder, "")
    s3_key = s3_folder + relative
    s3_client.upload_file(file, bucket_name, s3_key)
    if not quiet:
        print(f"Uploaded {file}")

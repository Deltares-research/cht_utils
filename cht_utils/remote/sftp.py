"""SSH/SFTP client for remote file operations."""

import os
import posixpath
import socket
from stat import S_ISDIR
from typing import Optional

import paramiko


class SSHSession:
    """SSH session with SFTP file transfer capabilities.

    Parameters
    ----------
    hostname : str
        Remote host address.
    username : str
        SSH username.
    key_file : str or None
        Path to private key file (not yet implemented).
    password : str or None
        SSH password.
    """

    def __init__(
        self,
        hostname: str,
        username: str = "root",
        key_file: Optional[str] = None,
        password: Optional[str] = None,
    ) -> None:
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.connect((hostname, 22))
        self.t = paramiko.Transport(self.sock)
        self.t.start_client()
        self.t.auth_password(username, password, fallback=False)
        self.sftp = paramiko.SFTPClient.from_transport(self.t)

    def command(self, cmd: str) -> str:
        """Execute a remote command.

        Parameters
        ----------
        cmd : str
            Command string (may contain newlines for multiple commands).

        Returns
        -------
        str
            Combined server response.
        """
        chan = self.t.open_session()
        chan.get_pty()
        chan.invoke_shell()
        chan.settimeout(20.0)
        ret = ""
        try:
            ret += chan.recv(1024).decode()
        except Exception:
            chan.send("\n")
            ret += chan.recv(1024).decode()
        for line in cmd.split("\n"):
            chan.send(line.strip() + "\n")
            ret += chan.recv(1024).decode()
        return ret

    def put(self, localfile: str, remotefile: str) -> None:
        """Upload a file.

        Parameters
        ----------
        localfile : str
            Local file path.
        remotefile : str
            Remote destination path.
        """
        self.sftp.put(localfile, remotefile)

    def put_all(self, localpath: str, remotepath: str) -> None:
        """Recursively upload a directory.

        Parameters
        ----------
        localpath : str
            Local directory to upload.
        remotepath : str
            Remote destination directory.
        """
        current_path = os.getcwd()
        os.chdir(os.path.split(localpath)[0])
        parent = os.path.split(localpath)[1]
        for walker in os.walk(parent):
            try:
                self.sftp.mkdir(remotepath + "/" + walker[0].replace("\\", "/"))
            except OSError:
                pass
            for file in walker[2]:
                remote_file = f"{remotepath}/{walker[0]}/{file}".replace("\\", "/")
                self.put(os.path.join(walker[0], file), remote_file)
        os.chdir(current_path)

    def get(self, remotefile: str, localfile: str) -> None:
        """Download a file.

        Parameters
        ----------
        remotefile : str
            Remote file path.
        localfile : str
            Local destination path.
        """
        self.sftp.get(remotefile, localfile)

    def sftp_walk(self, remotepath: str):
        """Walk a remote directory tree (like ``os.walk``).

        Parameters
        ----------
        remotepath : str
            Remote directory path.

        Yields
        ------
        tuple of (str, list, list)
            ``(path, folders, files)`` for each directory level.
        """
        files = []
        folders = []
        for f in self.sftp.listdir_attr(remotepath):
            if S_ISDIR(f.st_mode):
                folders.append(f.filename)
            else:
                files.append(f.filename)
        yield remotepath, folders, files
        for folder in folders:
            new_path = os.path.join(remotepath, folder).replace("\\", "/")
            yield from self.sftp_walk(new_path)

    def get_all(self, remotepath: str, localpath: str) -> None:
        """Recursively download a directory.

        Parameters
        ----------
        remotepath : str
            Remote directory path.
        localpath : str
            Local destination directory.
        """
        self.sftp.chdir(os.path.split(remotepath)[0])
        parent = os.path.split(remotepath)[1]
        os.makedirs(localpath, exist_ok=True)
        for walker in self.sftp_walk(parent):
            local_dir = os.path.join(localpath, walker[0])
            os.makedirs(local_dir, exist_ok=True)
            for file in walker[2]:
                self.get(
                    os.path.join(walker[0], file),
                    os.path.join(localpath, walker[0], file),
                )

    def write_command(self, text: str, remotefile: str) -> None:
        """Write text to a remote file and make it executable.

        Parameters
        ----------
        text : str
            Script content.
        remotefile : str
            Remote file path.
        """
        self.sftp.open(remotefile, "w").write(text)
        self.sftp.chmod(remotefile, 0o755)

    def rmtree(self, remotepath: str, level: int = 0) -> None:
        """Recursively delete a remote directory tree.

        Parameters
        ----------
        remotepath : str
            Remote directory to remove.
        level : int
            Recursion depth (internal use).
        """
        for f in self.sftp.listdir_attr(remotepath):
            rpath = posixpath.join(remotepath, f.filename)
            if S_ISDIR(f.st_mode):
                self.rmtree(rpath, level=level + 1)
            else:
                self.sftp.remove(rpath)
        self.sftp.rmdir(remotepath)

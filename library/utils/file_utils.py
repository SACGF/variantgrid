import gzip
import logging
import os
import subprocess
from collections.abc import Iterable
from dataclasses import dataclass
from hashlib import md5
from pathlib import Path
from typing import Optional, Union, IO


def open_file_or_filename(f, mode='r', **kwargs):
    """ If f is a string (filename), opens (handling gzip)
        otherwise does nothing if already a file object """
    if isinstance(f, str):
        if 'w' in mode:  # Create path if writing
            mk_path_for_file(f)

        return open_handle_gzip(f, mode, **kwargs)
    if all(hasattr(f, method) for method in ["read", "readlines"]):
        return f  # Already a File object
    raise ValueError(f"'{f}' ({type(f)}) not a file or string")


def mk_path(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def mk_path_for_file(f):
    mk_path(os.path.dirname(f))


def rm_if_exists(path):
    if os.path.exists(path):
        os.unlink(path)


def name_from_filename(filename: str, remove_gz=False) -> str:
    """Gets file name without extension or directory"""
    if remove_gz:
        filename = remove_gz_if_exists(filename)
    name = os.path.splitext(os.path.basename(filename))[0]
    return name


def file_to_array(filename, comment: Optional[str] = None, max_lines: Optional[int] = None):
    array = []
    f_or_f = open_file_or_filename(filename)
    for i, line in enumerate(f_or_f):
        if max_lines is not None and i >= max_lines:
            break

        if comment and line.startswith(comment):
            continue
        array.append(line.rstrip())
    return array


def file_or_filename_md5sum(file_or_filename: Union[IO, str]) -> str:
    # MD5 is used here only as a fast file-equivalency checksum (detecting identical
    # file contents), NOT for any security/cryptographic purpose, so collision
    # resistance is not required.
    m = md5()
    f: IO
    if hasattr(file_or_filename, "read"):
        f = file_or_filename  # Already a file
    else:
        f = open(file_or_filename, "rb")

    for chunk in iter(lambda: f.read(8192), b''):
        m.update(chunk)
    return m.hexdigest()


def remove_gz_if_exists(filename):
    GZIP_EXTENSIONS = [".gz", ".bgz"]
    for ext in GZIP_EXTENSIONS:
        if filename.endswith(ext):
            filename = os.path.splitext(filename)[0]  # remove extension
            break
    return filename


def get_extension_without_gzip(filename: str) -> str:
    filename = remove_gz_if_exists(filename)
    (_, ext) = os.path.splitext(filename)
    return ext[1:]


class IteratorFile:
    """ By Michael Anderson: http://stackoverflow.com/a/12593795 """

    def __init__(self, it):
        self.it = it
        self.next_chunk = ""

    def _grow_chunk(self):
        self.next_chunk = self.next_chunk + next(self.it)

    def read(self, n: int) -> str:
        if self.next_chunk is None:
            return ''
        try:
            while len(self.next_chunk) < n:
                self._grow_chunk()
            rv = self.next_chunk[:n]
            self.next_chunk = self.next_chunk[n:]
            return rv
        except StopIteration:
            rv = self.next_chunk
            self.next_chunk = None
            return rv

    def readline(self) -> Optional[str]:
        if self.next_chunk is None:
            return None
        try:
            while "\n" not in self.next_chunk:
                self._grow_chunk()
            n = self.next_chunk.index('\n') + 1
            rv = self.next_chunk[:n]
            self.next_chunk = self.next_chunk[n:]
            return rv
        except StopIteration:
            rv = self.next_chunk
            self.next_chunk = None
            return rv


class StashFile:
    """ File-like object that holds a value.

        Parts are accumulated and joined on read - repeated string concatenation goes quadratic
        once many writes stack up before a read (as when buffering rows for export) """

    def __init__(self):
        self.parts = []

    def write(self, value):
        self.parts.append(value)

    @property
    def value(self):
        data = "".join(self.parts)
        self.parts = []
        return data


def add_permissions_to_file(filename: str, add_stat: int):
    """ Adds file permission on a existing file path """
    st = os.stat(filename)
    try:
        os.chmod(filename, st.st_mode | add_stat)
    except Exception as e:
        logging.debug("Path '%s' stat is %s", filename, st)
        raise e


def open_handle_gzip(filename: str, mode=None, **kwargs):
    if filename.endswith(".gz"):
        open_func = gzip.open
    else:
        open_func = open

    args = ()
    if mode:
        args = (mode,)
    return open_func(filename, *args, **kwargs)


@dataclass
class DiskUsage:
    """ Free space on one mounted filesystem. Callers pass the minimum they need, so one class serves
        every threshold - SERVER_MIN_DISK_WARNING_GIGS for deployment warnings,
        ANNOTATION_MIN_FREE_DISK_GIGS for the annotation dispatcher, and so on. """
    mount_point: str
    available_kb: int
    available_nice: str  # human readable
    percent_nice: str  # human readable

    def has_capacity(self, min_gigs: float) -> bool:
        return self.available_kb >= min_gigs * 1000000

    def as_status_message(self, min_gigs: float) -> tuple[str, str]:
        """ (status, message) describing this mount against min_gigs. """
        message = f"Mount point '{self.mount_point}' ({self.percent_nice} used, {self.available_nice} available)"
        if self.has_capacity(min_gigs):
            return "info", message
        return "warning", message + f" is below minimum of {min_gigs}G"


def get_mount_point_for_directory(directory: str, mount_points: Iterable[str]) -> Optional[str]:
    """ The mount point `directory` lives on - the longest one that contains it.

        Longest wins because every absolute path is under '/', so a plain prefix test attributes
        '/data/annotation_scratch' to whichever of '/' and '/data' df happened to list first. Matching on
        a trailing separator also keeps '/database' from being read as living under '/data'.

        Purely textual, so it answers for a directory that does not exist yet. """
    best = None
    for mount_point in mount_points:
        prefix = mount_point if mount_point.endswith(os.sep) else mount_point + os.sep
        if directory == mount_point or directory.startswith(prefix):
            if best is None or len(mount_point) > len(best):
                best = mount_point
    return best


def get_disk_usage_for_directory(directory: str) -> Optional[DiskUsage]:
    """ Usage of the filesystem holding `directory`, or None if it can't be determined. """
    disk_usage = get_disk_usage()
    mount_point = get_mount_point_for_directory(directory, disk_usage)
    if mount_point is None:
        return None
    data = disk_usage[mount_point]
    nice_disk_usage = get_disk_usage(human_readable=True)
    return DiskUsage(mount_point=mount_point,
                     available_kb=int(data["avail"]),
                     available_nice=nice_disk_usage[mount_point]["avail"],
                     percent_nice=data["percent"])


def get_disk_usage(human_readable=False):
    """ returns a dict of DF output by mounted on """

    DF_COLUMNS = ["filesystem", "size", "used", "avail", "percent", "mounted_on"]
    cmd = ["df"]
    if human_readable:
        cmd.append("-h")

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdout_bytes, _ = p.communicate()
    stdout = stdout_bytes.decode()
    lines = stdout.split("\n")[1:]  # Skip header

    disk_percent_used = {}
    for line in lines:
        if line:
            columns = line.split()
            data = dict(zip(DF_COLUMNS, columns))
            mounted_on = data.pop("mounted_on")
            disk_percent_used[mounted_on] = data
    return disk_percent_used

import gzip
from hashlib import md5
import logging
import os
from pathlib import Path
import subprocess


def file_or_filename(f, mode='r'):
    if isinstance(f, str):
        if 'w' in mode:  # Create path if writing
            mk_path_for_file(f)

        return open(f, mode)
    if all([hasattr(f, method) for method in ["read", "readlines"]]):
        return f  # Already a File object
    raise ValueError(f"'{f}' ({type(f)}) not a file or string")


def mk_path(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def mk_path_for_file(f):
    mk_path(os.path.dirname(f))


def rm_if_exists(path):
    if os.path.exists(path):
        os.unlink(path)


def name_from_filename(filename, remove_gz=False):
    """Gets file name without extension or directory"""
    if remove_gz:
        filename = remove_gz_if_exists(filename)
    name = os.path.splitext(os.path.basename(filename))[0]
    return name


def file_to_array(filename, comment=None, max_lines=None):
    array = []
    f_or_f = file_or_filename(filename)
    for i, line in enumerate(f_or_f):
        if max_lines is not None and i > max_lines:
            break

        if comment and line.startswith(comment):
            continue
        array.append(line.rstrip())
    return array


def file_md5sum(filename):
    m = md5()
    with open(filename, "rb") as f:
        m.update(f.read())
    return m.hexdigest()


def remove_gz_if_exists(filename):
    GZIP_EXTENSIONS = [".gz", ".bgz"]
    for ext in GZIP_EXTENSIONS:
        if filename.endswith(ext):
            filename = os.path.splitext(filename)[0]  # remove extension
            break
    return filename


def get_extension_without_gzip(filename):
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

    def read(self, n):
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

    def readline(self):
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


def add_permissions_to_file(filename, add_stat):
    """ Adds file permission on a existing file path """
    st = os.stat(filename)
    try:
        os.chmod(filename, st.st_mode | add_stat)
    except Exception as e:
        logging.error("Path '%s' stat is %s", filename, st)
        raise e


def open_handle_gzip(filename, mode=None):
    if filename.endswith(".gz"):
        open_func = gzip.open
    else:
        open_func = open

    args = ()
    if mode:
        args = (mode,)
    return open_func(filename, *args)


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

import fcntl
import os
import stat

import nltk
from nltk.downloader import Downloader


def _remove_group_and_other_write(download_dir):
    """ nltk refuses to load data from a directory others can write to, and its downloader creates them under the
        process umask (often 002). Tightens whatever this account owns """
    uid = os.getuid()
    paths = [download_dir]
    for dirpath, dirnames, filenames in os.walk(download_dir):
        paths.extend(os.path.join(dirpath, name) for name in dirnames + filenames)
    for path in paths:
        st = os.lstat(path)
        if st.st_uid == uid and not stat.S_ISLNK(st.st_mode) and st.st_mode & (stat.S_IWGRP | stat.S_IWOTH):
            os.chmod(path, stat.S_IMODE(st.st_mode) & ~(stat.S_IWGRP | stat.S_IWOTH))


def ensure_nltk_data(resource):
    """ Downloads on first use. Locked as concurrent processes (celery workers, parallel tests) would otherwise
        read a corpus another is still writing - find() succeeds as soon as the partial zip exists """
    download_dir = Downloader().default_download_dir()
    os.makedirs(download_dir, exist_ok=True)
    with open(os.path.join(download_dir, ".download.lock"), "a", encoding="utf-8") as lock_file:
        fcntl.flock(lock_file, fcntl.LOCK_EX)
        _remove_group_and_other_write(download_dir)
        try:
            nltk.data.find(resource)
        except LookupError:
            nltk.download(resource.split('/')[-1], download_dir=download_dir)
            _remove_group_and_other_write(download_dir)

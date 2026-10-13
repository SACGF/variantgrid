import fcntl
import os

import nltk
from nltk.downloader import Downloader


def ensure_nltk_data(resource):
    """ Downloads on first use. Locked as concurrent processes (celery workers, parallel tests) would otherwise
        read a corpus another is still writing - find() succeeds as soon as the partial zip exists """
    download_dir = Downloader().default_download_dir()
    os.makedirs(download_dir, exist_ok=True)
    with open(os.path.join(download_dir, ".download.lock"), "a", encoding="utf-8") as lock_file:
        fcntl.flock(lock_file, fcntl.LOCK_EX)
        try:
            nltk.data.find(resource)
        except LookupError:
            nltk.download(resource.split('/')[-1], download_dir=download_dir)

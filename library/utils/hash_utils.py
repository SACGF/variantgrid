import hashlib
import json
import secrets
from hashlib import md5

import deprecation


def string_deterministic_hash(s: str) -> int:
    """
    When you want the same string to always hash to the same value (where security
    isn't a concern). Hashlib seemed overkill for this purpose.
    """
    val = 0
    for c in s:
        val = val * 13 + ord(c)
    return val


def _hash_str(method: callable, s: str) -> str:
    s_bytes = s.encode()
    return method(s_bytes).hexdigest()


@deprecation.deprecated(details="Use sha256sum_str instead")
def md5sum_str(s: str) -> str:
    return _hash_str(hashlib.md5, s)


def sha256sum_str(s: str) -> str:
    return _hash_str(hashlib.sha256, s)


def file_sha256sum(filename: str) -> str:
    hasher = hashlib.sha256()
    with open(filename, "rb") as f:
        hasher.update(f.read())
    return hasher.hexdigest()


def stable_dict_hash(d: dict) -> str:
    """ Order-independent sha256 of a JSON-serialisable dict, returned as a hex str (so it
        compares equal to a value round-tripped through a CharField). """
    s = json.dumps(d, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(s.encode()).hexdigest()


def secure_random_string() -> str:
    random_bytes = secrets.token_bytes(32)
    return hashlib.sha256(random_bytes).hexdigest()


@deprecation.deprecated(details="Use file_sha256sum instead")
def file_md5sum(filename: str):
    m = md5()
    with open(filename, "rb") as f:
        m.update(f.read())
    return m.hexdigest()

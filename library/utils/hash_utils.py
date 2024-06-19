import hashlib


def string_deterministic_hash(s: str) -> int:
    """
    When you want the same string to always hash to the same value (where security
    isnt a concern). Hashlib seemed overkill for this purpose.
    """
    val = 0
    for c in s:
        val = val * 13 + ord(c)
    return val


def _hash_str(method: callable, s: str) -> str:
    s_bytes = s.encode()
    return method(s_bytes).hexdigest()


def md5sum_str(s: str) -> str:
    return _hash_str(hashlib.md5, s)


def sha256_str(s: str) -> str:
    return _hash_str(hashlib.sha256, s)

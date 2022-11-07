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


def md5sum_str(s: str) -> str:
    s_bytes = s.encode()
    return hashlib.md5(s_bytes).hexdigest()


def sha1_str(s: str) -> str:
    s_bytes = s.encode()
    return hashlib.sha1(s_bytes).hexdigest()

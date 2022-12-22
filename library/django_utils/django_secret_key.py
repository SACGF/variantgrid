import os
import random
from base64 import urlsafe_b64encode as b64encode

random.seed()


def generate_key(max_length, seed_length):
    """
    Generate a Base64-encoded 'random' key by hashing the data.
    data is a tuple of seeding values. Pass arbitrary encoder and
    digester for specific hashing and formatting of keys

    From: https://gist.github.com/airtonix/6204802

    """
    PATTERN = "%%0%dX"
    JUNK_LEN = 1024
    junk = (PATTERN % (JUNK_LEN * 2)) % random.getrandbits(JUNK_LEN * seed_length)
    key = str(junk).encode()
    return b64encode(key)[:max_length]


def get_or_create_django_secret_key(key_dir):
    key_filename = os.path.join(key_dir, "django_secret_key.txt")
    if not os.path.exists(key_filename):
        secret_key = generate_key(50, 128)
        with open(key_filename, "wb") as f:
            f.write(secret_key)
    else:
        with open(key_filename, encoding="utf-8") as f:
            secret_key = f.read().strip()

    return secret_key

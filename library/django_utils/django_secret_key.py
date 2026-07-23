import os

from django.core.management.utils import get_random_secret_key


def get_or_create_django_secret_key(key_dir):
    key_filename = os.path.join(key_dir, "django_secret_key.txt")
    if not os.path.exists(key_filename):
        secret_key = get_random_secret_key()
        with open(key_filename, "w", encoding="utf-8") as f:
            f.write(secret_key)
    else:
        with open(key_filename, encoding="utf-8") as f:
            secret_key = f.read().strip()

    return secret_key

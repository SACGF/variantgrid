#!/usr/bin/env python3
import os
import sys

if __name__ == "__main__":
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "variantgrid.settings")

    from django.core.management import execute_from_command_line

    if sys.version_info < (3, 10):
        raise SystemExit("VariantGrid requires Python 3.10 or later.")

    execute_from_command_line(sys.argv)

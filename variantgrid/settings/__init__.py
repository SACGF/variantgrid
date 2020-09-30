"""
    This works similar to how Django Split settings works, except PyDev
    didn't properly handle split_settings.include and so had undefined variables
    (hence I import via direct code below)

    See https://code.djangoproject.com/wiki/SplitSettings for a discussion of techniques
"""

import logging
import os
import re
import socket

DEFAULT_DJANGO_SETTINGS = "variantgrid.settings"

django_settings_module = os.environ.get("DJANGO_SETTINGS_MODULE")
if django_settings_module == DEFAULT_DJANGO_SETTINGS:
    # Default settings, so do defaults then hostname-based .
    """
    (defaults are now loaded at the top of each setting file)
    from variantgrid.settings.components.default_settings import *
    from variantgrid.settings.components.celery_settings import *
    from variantgrid.settings.components.seqauto_settings import *
    """
    flattened_hostname = socket.gethostname().lower().split('.')[0].replace('-', '')
    pwd = os.path.dirname(__file__)
    flattened_hostname_path = os.path.join(pwd, 'env', f"{flattened_hostname}.py")
    logging.info(f'LOADING settings file {flattened_hostname_path}')
    if os.path.exists(flattened_hostname_path):
        flattened_hostname_module = f"variantgrid.settings.env.{flattened_hostname}"
        exec(f"from {flattened_hostname_module} import *")
else:
    exec(f"from {django_settings_module} import *")

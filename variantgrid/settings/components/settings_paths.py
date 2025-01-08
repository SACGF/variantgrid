# The rest of the settings files imports this to use as base of paths

import os

# Up 2 more dirs than normal (as we're in variantgrid.settings.components dir)
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SETTINGS_DIR = os.path.dirname(_THIS_DIR)
BASE_DIR = os.path.dirname(os.path.dirname(SETTINGS_DIR))

PRIVATE_DATA_ROOT = os.path.join(BASE_DIR, "data")

VARIANTGRID_REPO_APP_DIR = os.path.join(BASE_DIR, "variantgrid")
VARIANTGRID_REPO_REFERENCE_DIR = os.path.join(VARIANTGRID_REPO_APP_DIR, "data", "reference")

ANNOTATION_BASE_DIR = "/data/annotation"


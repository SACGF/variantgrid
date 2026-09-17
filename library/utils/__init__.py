"""
The `library.utils` facade: every library/utils/*_utils.py module is star-imported here so callers write
`from library.utils import first, batch_iterator, sha256sum_str`. Check the per-topic modules
(collection, text, html, json, hash, file, date, diff, export, model, class, os, timer, misc, color)
before adding a helper - library/CLAUDE.md lists what each holds. A new utils module must be added
to this list to be reachable through the facade.
"""
from library.utils.class_utils import *
from library.utils.collection_utils import *
from library.utils.color_utils import *
from library.utils.date_utils import *
from library.utils.diff_utils import *
from library.utils.export_utils import *
from library.utils.file_utils import *
from library.utils.hash_utils import *
from library.utils.html_utils import *
from library.utils.json_utils import *
from library.utils.misc_utils import *
from library.utils.model_utils import *
from library.utils.os_utils import *
from library.utils.text_utils import *
from library.utils.timer_utils import *

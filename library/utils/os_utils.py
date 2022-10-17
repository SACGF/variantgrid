import logging
import subprocess
from typing import Tuple, Optional

from django.conf import settings


def execute_cmd(cmd: list, **kwargs) -> Tuple[int, Optional[str], Optional[str]]:
    shell = kwargs.get("shell", settings.POPEN_SHELL)

    if shell:
        command = ' '.join(cmd)
        logging.info('About to call %s', command)
        pipes = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)
        logging.info('Completed')
    else:
        pipes = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)

    std_out, std_err = pipes.communicate()
    return pipes.returncode, std_out.decode() if std_out else None, std_err.decode() if std_err else None


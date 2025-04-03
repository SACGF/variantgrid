import logging
import subprocess
from dataclasses import dataclass
from typing import Optional

from library.utils import FormerTuple


@dataclass(frozen=True)
class CmdOutput(FormerTuple):
    return_code: int
    std_out: Optional[str]
    std_err: Optional[str]

    @property
    def as_tuple(self) -> tuple:
        return self.return_code, self.std_out, self.std_err


def execute_cmd(cmd: list, **kwargs) -> CmdOutput:
    if kwargs.pop("shell", False):
        command = ' '.join(cmd)
        logging.info('About to call %s', command)
        pipes = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)
        logging.info('Completed')
    else:
        logging.info('About to call %s', cmd)
        pipes = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)

    std_out, std_err = pipes.communicate()
    return CmdOutput(
        return_code=pipes.returncode,
        std_out=std_out.decode() if std_out else None,
        std_err=std_err.decode() if std_err else None
    )

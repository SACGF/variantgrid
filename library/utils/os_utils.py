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
    # The shell=True branch is opt-in and currently used by exactly one caller
    # (snpdb.bcftools_liftover) to build a trusted, server-side command pipeline -
    # no user-supplied input reaches `cmd`. Keep it that way: never pass
    # user-controlled values with shell=True.
    #
    # process_callback (#1658): optional callable invoked with the live Popen right after it starts,
    # so a caller (e.g. the annotation lease heartbeat on a separate thread) can .kill() the subprocess
    # while communicate() blocks here - the resulting non-zero return_code lets the caller abort cleanly.
    process_callback = kwargs.pop("process_callback", None)
    if kwargs.pop("shell", False):
        command = ' '.join(cmd)
        logging.info('About to call %s', command)
        pipes = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)
        logging.info('Completed')
    else:
        logging.info('About to call %s', cmd)
        pipes = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)

    if process_callback is not None:
        process_callback(pipes)

    std_out, std_err = pipes.communicate()
    return CmdOutput(
        return_code=pipes.returncode,
        std_out=std_out.decode() if std_out else None,
        std_err=std_err.decode() if std_err else None
    )

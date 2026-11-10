import logging
import os
import signal
import subprocess
import threading
import time
from dataclasses import dataclass
from typing import Optional

from library.utils.collection_utils import FormerTuple


@dataclass(frozen=True)
class CmdOutput(FormerTuple):
    return_code: int
    std_out: Optional[str]
    std_err: Optional[str]
    # Peak resident memory of the largest single process the command ran (itself or a descendant it
    # waited for). Not part of the tuple, so `return_code, std_out, std_err = execute_cmd(...)` still works
    max_rss_mb: Optional[int] = None

    @property
    def as_tuple(self) -> tuple:
        return self.return_code, self.std_out, self.std_err


def format_called_process_error(e: subprocess.CalledProcessError) -> str:
    """ str(CalledProcessError) is only the command and exit status - tools like bcftools put the
        diagnostics the user actually needs on stderr, so include whatever was captured """

    def _decode(output) -> str:
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        return (output or "").strip()

    message = str(e)
    for label, output in (("stderr", _decode(e.stderr)), ("stdout", _decode(e.stdout))):
        if output:
            message += f"\n{label}: {output}"
    return message


def kill_process_tree(process: subprocess.Popen):
    """ SIGKILL a process started by execute_cmd along with everything it started: execute_cmd gives each
        command its own process group, so a tool's children (AnnotSV's bedtools, VEP's forks, the later
        stages of a shell pipeline) go with it rather than outliving the run that started them (#2749) """
    if process.returncode is None:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass


def _communicate_with_rusage(process: subprocess.Popen, timeout: Optional[float]):
    """ Popen.communicate, except the child is reaped with os.wait4 - the only call that hands back its
        resource usage, which Popen does not keep. Returns (std_out, std_err, rusage, timed_out); rusage is
        None if something else reaped the child first. """
    outputs = {}

    def _read(name, pipe):
        outputs[name] = pipe.read()
        pipe.close()

    readers = [threading.Thread(target=_read, args=(name, pipe), daemon=True)
               for name, pipe in (("std_out", process.stdout), ("std_err", process.stderr))]
    for reader in readers:
        reader.start()

    deadline = None if timeout is None else time.monotonic() + timeout
    timed_out = False
    rusage = None
    try:
        for reader in readers:
            reader.join(None if deadline is None else max(0.0, deadline - time.monotonic()))
        wait_flags = 0 if deadline is None else os.WNOHANG
        while True:
            pid, wait_status, rusage = os.wait4(process.pid, wait_flags)
            if pid == process.pid:
                process.returncode = os.waitstatus_to_exitcode(wait_status)
                break
            if time.monotonic() >= deadline:
                # Kill then drain - the pipes only reach EOF once the whole tree has let go of them
                timed_out = True
                kill_process_tree(process)
                wait_flags = 0
            else:
                time.sleep(0.05)
    except ChildProcessError:
        rusage = None
        process.wait()
    for reader in readers:
        reader.join()
    return outputs["std_out"], outputs["std_err"], rusage, timed_out


def execute_cmd(cmd: list, **kwargs) -> CmdOutput:
    # The shell=True branch is opt-in and currently used by exactly one caller
    # (snpdb.bcftools_liftover) to build a trusted, server-side command pipeline -
    # no user-supplied input reaches `cmd`. Keep it that way: never pass
    # user-controlled values with shell=True.
    #
    # process_callback (#1658): optional callable invoked with the live Popen right after it starts,
    # so a caller (e.g. the annotation lease heartbeat on a separate thread) can kill_process_tree() it
    # while this blocks - the resulting non-zero return_code lets the caller abort cleanly.
    process_callback = kwargs.pop("process_callback", None)
    # timeout (#720): seconds to wait before killing the subprocess, as subprocess.run's timeout. Callers
    # get whatever was written before the kill plus a non-zero return_code, so a timeout needs no branch
    # of its own on top of the return_code handling they already have. Default None = wait forever.
    timeout = kwargs.pop("timeout", None)
    # Its own session, so its own process group: what kill_process_tree signals
    kwargs["start_new_session"] = True
    if kwargs.pop("shell", False):
        command = ' '.join(cmd)
        logging.info('About to call %s', command)
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)
    else:
        logging.info('About to call %s', cmd)
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)

    if process_callback is not None:
        process_callback(process)

    std_out, std_err, rusage, timed_out = _communicate_with_rusage(process, timeout)
    if timed_out:
        logging.warning("Command timed out after %ss: %s", timeout, cmd)

    return CmdOutput(
        return_code=process.returncode,
        std_out=std_out.decode() if std_out else None,
        std_err=std_err.decode() if std_err else None,
        max_rss_mb=rusage.ru_maxrss // 1024 if rusage else None,  # Linux reports KB
    )

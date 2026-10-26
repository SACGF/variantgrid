"""
Whether the installed systemd units can start on this checkout. The units in `config/systemd/` are the same
file for every install location: they read VG_INSTALL_DIR from /etc/variantgrid/variantgrid.env and the unit's
own .env (`config/systemd/README.txt`). A deployment that copied a unit without those files, or whose
VG_INSTALL_DIR points at another checkout, fails at the next restart - so this runs before it, in
`manage.py deployment_check`, and stops the upgrade. So does a unit's own .env left over from an older unit:
the env files are copied with `cp -rn`, keeping a host's tuning, so one that no longer suits its unit stays.

Entry points: `check_systemd_units`, run by `deployment_check`; `systemd_env_files_problem` is the predicate the
manual migration that installs the env files shares; `detached_simple_units` finds a leftover .env a Type=simple unit
can't run with; `systemctl_show` is the one place a unit's properties are read.
"""
import os
import shlex
import subprocess

from django.conf import settings

from library.vg.services import service_names

ETC_DIR = "/etc/variantgrid"
SHARED_ENV_FILE = f"{ETC_DIR}/variantgrid.env"
INSTALL_STEP = ("install the env files: sudo cp -rn config/celery config/gunicorn config/variantgrid.env "
                f"{ETC_DIR}/ then set VG_INSTALL_DIR in {SHARED_ENV_FILE}, "
                "sudo cp config/systemd/*.service /lib/systemd/system/ and sudo systemctl daemon-reload")


def systemctl_show(unit: str, *properties: str) -> dict[str, list[str]] | None:
    """ The unit's properties, each as the list of lines systemctl printed for it (EnvironmentFiles is one line
        per file), or None when there is no such loaded unit: a dev box running runserver, a host without systemd """
    args = ["systemctl", "show", unit, "-p", "LoadState"]
    for prop in properties:
        args += ["-p", prop]
    try:
        result = subprocess.run(args, capture_output=True, text=True, check=False)
    except FileNotFoundError:
        return None
    if result.returncode != 0:
        return None
    props: dict[str, list[str]] = {}
    for line in result.stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            props.setdefault(key, []).append(value)
    if props.get("LoadState") != ["loaded"]:
        return None
    return props


def _environment_file_paths(props: dict[str, list[str]]) -> list[str]:
    """ Paths of the unit's EnvironmentFile= lines, from "/path (ignore_errors=yes)" as systemctl prints them """
    return [value.split(" (", 1)[0] for value in props.get("EnvironmentFiles", []) if value]


def unit_environment_files(unit: str) -> list[str] | None:
    """ Paths of the unit's EnvironmentFile= lines, None when the unit isn't installed """
    props = systemctl_show(unit, "EnvironmentFiles")
    if props is None:
        return None
    return _environment_file_paths(props)


def env_file_value(path: str, name: str) -> str | None:
    """ `name=` from an env file with its quotes removed, as systemd reads it; None when the file or the
        variable is missing """
    try:
        with open(path, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line.startswith(f"{name}="):
                    value = shlex.split(line.split("=", 1)[1])
                    return value[0] if value else None
    except FileNotFoundError:
        return None
    return None


def shared_env_install_dir() -> str | None:
    """ VG_INSTALL_DIR from the shared env file, None when the file or the variable is missing """
    return env_file_value(SHARED_ENV_FILE, "VG_INSTALL_DIR")


def detached_simple_units() -> dict[str, str]:
    """ {unit: env file} for each Type=simple unit whose own .env passes --detach in CELERYD_OPTS: the process
        systemd started exits once celery has forked, systemd stops the unit, kills the fork and restarts it
        until start-limit-hit. celeryd_beat became Type=simple without --detach, and an older .env survives
        the `cp -rn` install """
    detached = {}
    for unit in service_names():
        props = systemctl_show(unit, "Type", "EnvironmentFiles")
        if props is None or props.get("Type") != ["simple"]:
            continue
        for path in _environment_file_paths(props):
            if path != SHARED_ENV_FILE and "--detach" in (env_file_value(path, "CELERYD_OPTS") or "").split():
                detached[unit] = path
    return detached


def _installed_units_env_files() -> dict[str, list[str]]:
    return {unit: files for unit in service_names() if (files := unit_environment_files(unit)) is not None}


def systemd_env_files_problem() -> str | None:
    """ Why the installed env-driven units would fail to start, or None: no unit installed, or every env file
        they read exists and VG_INSTALL_DIR is this checkout. Units predating variantgrid.env (hardcoded paths)
        read no env file we know of, so they are not a problem here - `check_systemd_units` warns about them """
    units = _installed_units_env_files()
    missing = sorted({path for files in units.values() for path in files if not os.path.exists(path)})
    if missing:
        return f"missing {', '.join(missing)}"
    if any(SHARED_ENV_FILE in files for files in units.values()):
        install_dir = shared_env_install_dir()
        if install_dir is None:
            return f"{SHARED_ENV_FILE} does not set VG_INSTALL_DIR"
        if os.path.realpath(install_dir) != os.path.realpath(settings.BASE_DIR):
            return f"VG_INSTALL_DIR={install_dir} is not this checkout ({settings.BASE_DIR})"
    return None


def check_systemd_units() -> dict:
    units = _installed_units_env_files()
    if not units:
        return {}

    env_check = {"valid": True}
    if problem := systemd_env_files_problem():
        env_check = {"valid": False, "fix": f"{problem} - the services will not restart, {INSTALL_STEP}"}

    # A unit still hardcoding its install dir works, so a warning: the env-driven unit is one file for every box
    hardcoded_check = {"valid": True}
    if hardcoded := sorted(unit for unit, files in units.items() if SHARED_ENV_FILE not in files):
        hardcoded_check["warning"] = (f"{', '.join(hardcoded)} do not read {SHARED_ENV_FILE} (units from before "
                                      f"config/variantgrid.env) - {INSTALL_STEP}")

    detach_check = {"valid": True}
    if detached := detached_simple_units():
        copies = "; ".join(f"sudo cp config/celery/{os.path.basename(path)} {path}" for path in detached.values())
        detach_check = {"valid": False,
                        "fix": f"{', '.join(detached)} will not stay running: their .env passes --detach, which a "
                               f"Type=simple unit cannot run with. Refresh it from the checkout ({copies}, "
                               f"re-applying any local tuning) then sudo systemctl restart {' '.join(detached)}"}

    return {
        "env files": env_check,
        "env-driven units": hardcoded_check,
        "env files match units": detach_check,
    }

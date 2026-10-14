"""
Whether the installed systemd units can start on this checkout. The units in `config/systemd/` are the same
file for every install location: they read VG_INSTALL_DIR from /etc/variantgrid/variantgrid.env and the unit's
own .env (`config/systemd/README.txt`). A deployment that copied a unit without those files, or whose
VG_INSTALL_DIR points at another checkout, fails at the next restart - so this runs before it, in
`manage.py deployment_check`, and stops the upgrade.

Entry points: `check_systemd_units`, run by `deployment_check`; `systemd_env_files_problem` is the predicate the
manual migration that installs the env files shares; `systemctl_show` is the one place a unit's properties are read.
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


def unit_environment_files(unit: str) -> list[str] | None:
    """ Paths of the unit's EnvironmentFile= lines ("/path (ignore_errors=yes)" as systemctl prints them),
        None when the unit isn't installed """
    props = systemctl_show(unit, "EnvironmentFiles")
    if props is None:
        return None
    return [value.split(" (", 1)[0] for value in props.get("EnvironmentFiles", []) if value]


def shared_env_install_dir() -> str | None:
    """ VG_INSTALL_DIR from the shared env file, None when the file or the variable is missing """
    try:
        with open(SHARED_ENV_FILE, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line.startswith("VG_INSTALL_DIR="):
                    value = shlex.split(line.split("=", 1)[1])
                    return value[0] if value else None
    except FileNotFoundError:
        return None
    return None


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

    return {
        "env files": env_check,
        "env-driven units": hardcoded_check,
    }

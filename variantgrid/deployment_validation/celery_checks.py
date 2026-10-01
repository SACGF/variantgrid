import importlib

from celery import Task
from django.conf import settings

from library.utils import import_class
from variantgrid.celery import app


def _check_celery_tasks(setting_name, task_names, is_task_class=False) -> dict:
    bad_task_routes = {}
    for class_or_module_path in task_names:
        try:
            if is_task_class:
                klass = import_class(class_or_module_path)
                if not isinstance(klass, Task):
                    bad_task_routes[class_or_module_path] = "Not a celery task"
            else:
                importlib.import_module(class_or_module_path)
        except Exception:
            bad_task_routes[class_or_module_path] = "Not found"

    if bad_task_routes:
        data = {
            "valid": False,
            "fix": f"Edit settings.{setting_name}: {bad_task_routes}",
        }
    else:
        data = {"valid": True}
    return data


def check_celery_tasks() -> dict:
    celery_tasks = {}

    celery_settings_are_tasks = {
        "CELERY_TASK_ROUTES": True,
        "CELERY_IMPORTS": False,
    }

    for celery_setting, is_task_class in celery_settings_are_tasks.items():
        celery_tasks[celery_setting] = _check_celery_tasks(celery_setting, getattr(settings, celery_setting),
                                                           is_task_class=is_task_class)

    return celery_tasks


def _worker_task_modules() -> set[str]:
    """ The modules a worker imports on its own: CELERY_IMPORTS plus autodiscovered <app>.tasks """
    return set(settings.CELERY_IMPORTS) | {f"{app_name}.tasks" for app_name in settings.INSTALLED_APPS}


def check_beat_schedule_tasks(beat_schedule=None) -> dict:
    """ Beat only sends a task name, so the worker must register it by importing the defining module itself.
        A module that is only reached through a urls/views import chain is registered by the system checks Celery
        runs at startup, which CELERY_SKIP_CHECKS=1 turns off - so we require the explicit route. """
    if beat_schedule is None:
        beat_schedule = app.conf.beat_schedule
    worker_task_modules = _worker_task_modules()
    bad_tasks = {}
    for schedule_name, entry in beat_schedule.items():
        task_name = entry["task"]
        module_name = task_name.rsplit(".", 1)[0]
        try:
            if not isinstance(import_class(task_name), Task):
                bad_tasks[schedule_name] = f"{task_name} is not a celery task"
                continue
        except Exception:
            bad_tasks[schedule_name] = f"{task_name} not found"
            continue
        if module_name not in worker_task_modules:
            bad_tasks[schedule_name] = f"{module_name} is not in settings.CELERY_IMPORTS"

    if bad_tasks:
        return {"valid": False, "fix": f"Beat schedule tasks a worker would not register: {bad_tasks}"}
    return {"valid": True}

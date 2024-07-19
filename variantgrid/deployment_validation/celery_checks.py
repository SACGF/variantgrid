import importlib

from celery import Task
from django.conf import settings

from library.utils import import_class


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
        except:
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

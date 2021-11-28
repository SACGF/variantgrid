import os

from django.conf import settings
from django.template.loader import render_to_string

PBS_TEMPLATE = os.path.join(settings.SEQAUTO_DIR, "templates", "pbs_script.template")


def get_dependency_flags(*args):
    deps = ':'.join(["$%s" % b.get_variable_name() for b in args if b])
    if deps:
        dependency_flags = f" -W depend=afterok:{deps} "
    else:
        dependency_flags = ''

    return dependency_flags


def create_pbs_script(name, filename, out_file, command, cores, mem, job_script_pk):
    pythonpaths = [settings.BASE_DIR]
    pbs_pythonpath = settings.SEQAUTO_SCRIPT_PARAMS.get("pythonpath")
    if pbs_pythonpath:
        pythonpaths.append(pbs_pythonpath)
    pythonpath = ':'.join(pythonpaths)

    paths = [os.path.join(settings.BASE_DIR, "scripts"),
             os.path.join(settings.SEQAUTO_DIR, "scripts")]
    path = ':'.join(paths)

    job_complete = settings.SEQAUTO_JOB_COMPLETE
    if job_complete and settings.SEQAUTO_VIRTUALENV_RUNNER:
        job_complete = f"{settings.SEQAUTO_VIRTUALENV_RUNNER} {job_complete}"

    context = {"name": name,
               "email": settings.SEQAUTO_PBS_EMAIL,
               "out": out_file,
               "mem": mem,
               "cores": cores,
               "queue": settings.SEQAUTO_PBS_QUEUE,
               "path": path,
               "pythonpath": pythonpath,
               "command": command,
               "job_script_complete": job_complete,
               "job_script_pk": job_script_pk}
    script_contents = render_to_string(PBS_TEMPLATE, context)

    with open(filename, 'w') as f:
        f.write(script_contents)

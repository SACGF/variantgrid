
from django.core.exceptions import PermissionDenied
from django.db.utils import IntegrityError
from django.http.response import (
    HttpResponse,
    HttpResponseServerError,
)
from django.shortcuts import render
from django.views.decorators.http import require_POST
from guardian.shortcuts import get_objects_for_user

from library.django_utils import (
    add_save_message,
)
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsMixin
from library.guardian_utils import DjangoPermission
from library.utils import import_class
from snpdb import forms
from snpdb.models import (
    VCF,
    Sample,
)
from snpdb.tasks.soft_delete_tasks import soft_delete_samples, soft_delete_vcfs


def _import_permission_class(class_name):
    """ class_name comes from the URL - resolve it to a class, tolerating anything that doesn't. """
    try:
        klass = import_class(class_name)
    except (ImportError, AttributeError, ValueError):
        klass = None
    return klass if isinstance(klass, type) else None


def _require_guardian_permission_class(class_name):
    """ For the sharing views (group_permissions / bulk): only models that implement Guardian
        object-level permissions can have their permissions viewed/edited. Without this an attacker
        could pass an arbitrary class_name and act on models whose can_write() defaults to True
        (e.g. collaborative Wikis). """
    klass = _import_permission_class(class_name)
    if not (klass and issubclass(klass, GuardianPermissionsMixin)):
        raise PermissionDenied(f"'{class_name}' is not a permission-managed class")
    return klass


def _get_writable_object(user, klass, primary_key):
    name = klass.__name__
    obj = klass.objects.get(pk=primary_key)

    if not obj.can_write(user):
        write_perm = DjangoPermission.perm(obj, DjangoPermission.WRITE)
        msg = f"You do not have permission {write_perm} needed to modify {name}"
        raise PermissionDenied(msg)

    return obj, name


def get_writable_class_object(user, class_name, primary_key):
    klass = _require_guardian_permission_class(class_name)
    return _get_writable_object(user, klass, primary_key)


def get_writable_class_objects(user, class_name):
    klass = _require_guardian_permission_class(class_name)
    name = klass.__name__
    write_perm = DjangoPermission.perm(klass, DjangoPermission.WRITE)
    qs = get_objects_for_user(user, write_perm, klass=klass, accept_global_perms=False)
    return qs, name


def group_permissions(request, class_name, primary_key):
    obj, name = get_writable_class_object(request.user, class_name, primary_key)

    try:
        # If object has "get_permission_object" it can delegate it.
        permission_obj = obj.get_permission_object()
        perm_obj_name = permission_obj.__class__.__name__
    except AttributeError:
        # Default is use itself
        permission_obj = obj
        perm_obj_name = name

    permission_forms = get_group_permission_forms(request, permission_obj)

    if request.method == 'POST':
        valid = all([pf.is_valid() for pf in permission_forms])
        if valid:
            for pf in permission_forms:
                pf.save()
        add_save_message(request, valid, f"{perm_obj_name} group permissions")

    get_listing_url = getattr(obj, "get_listing_url", None)
    if get_listing_url:
        delete_redirect_url = get_listing_url()
    else:
        delete_redirect_url = "/"

    context = {'permission_forms': permission_forms,
               'class_name': class_name,
               'name': name,
               'perm_obj_name': perm_obj_name,
               'permission_obj': permission_obj,
               'instance': obj,
               'delete_redirect_url': delete_redirect_url}
    return render(request, 'snpdb/data/group_permissions.html', context)


SOFT_DELETE_FUNCS = {
    VCF: soft_delete_vcfs,
    Sample: soft_delete_samples,
}


@require_POST
def group_permissions_object_delete(request, class_name, primary_key):
    klass = _import_permission_class(class_name)
    if soft_delete := SOFT_DELETE_FUNCS.get(klass):
        # These cascade into a lot of data, so are marked for deletion then removed by a Celery task
        soft_delete(request.user, primary_key)
    else:
        # Deletion via this generic endpoint is opt-in per class - having WRITE permission isn't
        # enough (e.g. ClassificationModification audit records must never be deletable here).
        allow_delete = getattr(klass, "allow_group_permission_delete", None)
        if not (callable(allow_delete) and allow_delete()):
            raise PermissionDenied(f"'{class_name}' does not allow deletion via group permissions")
        obj, _ = _get_writable_object(request.user, klass, primary_key)
        try:
            # Some objects can't be removed once they've been used (eg AnalysisTemplate) so soft delete instead
            delete = getattr(obj, "delete_or_soft_delete", obj.delete)
            delete()
        except IntegrityError as ie:
            pks = ", ".join(str(o.pk) for o in ie.args[1])
            error_message = f"{ie.args[0]}: {pks}"
            return HttpResponseServerError(content=error_message)

    return HttpResponse()


def bulk_group_permissions(request, class_name):
    qs, name = get_writable_class_objects(request.user, class_name)
    groups = list(request.user.groups.all().order_by("name"))

    objects_and_forms = []
    for obj in qs:
        permission_forms = get_group_permission_forms(request, obj, groups=groups)
        objects_and_forms.append((obj, permission_forms))

    if request.method == 'POST':
        all_forms = []
        for _, permission_forms in objects_and_forms:
            all_forms.extend(permission_forms)
        valid = all([pf.is_valid() for pf in all_forms])
        if valid:
            for pf in all_forms:
                pf.save()
        add_save_message(request, valid, f"{name} group permissions")

    context = {"name": name,
               "groups": groups,
               "objects_and_forms": objects_and_forms}
    return render(request, 'snpdb/data/bulk_group_permissions.html', context)


def get_group_permission_forms(request, obj, groups=None):
    if groups is None:
        groups = request.user.groups.all().order_by("name")
    return [forms.GroupPermissionForm(request.POST or None, obj=obj, group=group) for group in groups]

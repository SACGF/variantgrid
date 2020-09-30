import csv
from typing import Optional

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
import json
import os

from snpdb.models.models import Lab, Organization
from snpdb.models.models_user_settings import UserSettings

module_dir = os.path.dirname(__file__)  # get current directory
lab_path = os.path.join(module_dir, 'labs.csv')
user_path = os.path.join(module_dir, 'users.csv')


def csv_to_dict(file_path):
    header = None
    with open(file_path, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        rows = []
        for csv_row in reader:
            if not header:
                header = csv_row
                print('header = ')
                print(header)
            else:
                row = {}
                for index, cell in enumerate(csv_row):
                    cell = cell.strip()
                    if cell.lower() == 'true':
                        cell = True
                    elif cell.lower() == 'false':
                        cell = False
                    elif len(cell) == 0:
                        cell = None
                    row[header[index]] = cell
                rows.append(row)
        return rows


def ensure_labs():
    lab_datas = csv_to_dict(lab_path)
    for lab_data in lab_datas:
        org_id = lab_data.pop('organization')
        lab_data['organization'] = Organization.objects.get(group_name=org_id)
        config = lab_data.pop('variant_classification_config')
        if config:
            config = json.loads(config)
            lab_data['variant_classification_config'] = config

        lab, created = Lab.objects.get_or_create(group_name=lab_data['group_name'], defaults=lab_data)
        if not created:
            Lab.objects.filter(pk=lab.id).update(**lab_data)


def ensure_users():
    user_datas = csv_to_dict(user_path)

    all_lab_groups = set()
    for lab in Lab.objects.all():
        all_lab_groups.add(lab.group)
        all_lab_groups.add(lab.group_institution)

    for user_data in user_datas:
        user, _ = User.objects.get_or_create(username=user_data['username'], defaults={
            'email': user_data['email'] or '',
            'first_name': user_data['first_name'] or '',
            'last_name': user_data['last_name'] or ''
        })
        lab_id: Optional[str] = user_data['lab']
        if lab_id is None:
            user.active = True
        elif lab_id == 'inactive':
            user.active = False
        else:
            user.active = True
            lab = Lab.objects.get(group_name=lab_id)
            for group in all_lab_groups:
                user.groups.remove(group)
            user.groups.add(lab.group)
            user.groups.add(lab.group_institution)
            # also add legacy labs
            parts = lab_id.split('/')
            org_part = parts[0]
            lab_part = parts[1]

            legacy_lab = Lab.objects.filter(group_name=org_part + '/legacy_' + lab_part).first()
            if legacy_lab:
                user.groups.add(legacy_lab.group)

            user_settings = UserSettings.get_for_user(user)
            if not user_settings.default_lab:
                user_settings.default_lab = lab
                user_settings.save()

        user.save()


class Command(BaseCommand):

    def handle(self, *args, **options):
        ensure_labs()
        ensure_users()

from typing import Tuple

from django.contrib.auth.models import User

from classification.models import Classification
from snpdb.models import Lab, Organization


class ClassificationTestUtils:

    @staticmethod
    def setUp():
        org = Organization(
            name='InstX',
            group_name='instx'
        )
        org.save()
        lab = Lab.objects.create(
            name='Labby',
            organization=org,
            city='CityA',
            country='CountryA',
            group_name='instx/labby'
        )

        ext_lab = Lab.objects.create(
            external=True,
            name='External',
            organization=org,
            city='CityB',
            country='CountryB',
            group_name='instx/ext'
        )

        user = User(
            username='joejoe',
            email='joe@joe.com'
        )
        user.save()
        user.groups.add(lab.group)
        user.save()

        user = User(
            username='joejoe2',
            email='joe@joe2.com'
        )
        user.save()
        user.groups.add(lab.group)
        user.groups.add(ext_lab.group)
        user.save()

    @staticmethod
    def tearDown():
        Lab.objects.filter(group_name='instx/labby').delete()
        Lab.objects.filter(group_name='instx/ext').delete()
        User.objects.filter(username__in=['joejoe', 'joejoe2']).delete()
        Classification.objects.all().delete()

    @staticmethod
    def lab_and_user() -> Tuple[Lab, User]:
        return Lab.objects.filter(group_name='instx/labby').get(), User.objects.filter(username='joejoe').get()

    @staticmethod
    def external_lab_and_user() -> Tuple[Lab, User]:
        return Lab.objects.filter(group_name='instx/ext').get(), User.objects.filter(username='joejoe').get()

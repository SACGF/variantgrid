from typing import Tuple

from django.contrib.auth.models import User

from snpdb.models import Lab, Organization
from classification.models import Classification


class ClassificationTestUtils:

    @staticmethod
    def setUp():
        org = Organization(
            name='InstX',
            group_name='instx'
        )
        org.save()
        lab = Lab(
            name='Labby',
            organization=org,
            city='CityA',
            country='CountryA',
            group_name='instx/labby'
        )
        lab.save()

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
        user.save()

    @staticmethod
    def tearDown():
        Lab.objects.filter(group_name='instx/labby').delete()
        User.objects.filter(username__in=['joejoe', 'joejoe2']).delete()
        Classification.objects.all().delete()

    @staticmethod
    def lab_and_user() -> Tuple[Lab, User]:
        return Lab.objects.filter(group_name='instx/labby').get(), User.objects.filter(username='joejoe').get()

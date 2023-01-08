from typing import Tuple

from django.contrib.auth.models import User

from classification.models import Classification
from snpdb.models import Lab, Organization, Country


class ClassificationTestUtils:

    @staticmethod
    def setUp():
        org = Organization.objects.create(
            name='InstX',
            group_name='instx'
        )
        country_a = Country.objects.get_or_create(name='CountryA')[0]
        country_b = Country.objects.get_or_create(name='CountryB')[0]

        lab = Lab.objects.create(
            name='Labby',
            organization=org,
            city='CityA',
            country=country_a,
            group_name='instx/labby'
        )

        ext_lab = Lab.objects.create(
            external=True,
            name='External',
            organization=org,
            city='CityB',
            country=country_b,
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
        pass
        # Lab.objects.filter(group_name='instx/labby').delete()
        # Lab.objects.filter(group_name='instx/ext').delete()
        # User.objects.filter(username__in=['joejoe', 'joejoe2']).delete()
        # Classification.objects.all().delete()

    @staticmethod
    def lab_and_user() -> Tuple[Lab, User]:
        return Lab.objects.get(group_name='instx/labby'), User.objects.get(username='joejoe')

    @staticmethod
    def external_lab_and_user() -> Tuple[Lab, User]:
        return Lab.objects.get(group_name='instx/ext'), User.objects.get(username='joejoe')

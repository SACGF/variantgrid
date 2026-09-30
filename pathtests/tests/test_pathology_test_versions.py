from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.urls import reverse
from django.utils import timezone

from genes.models import GeneList, GeneListGeneSymbol, GeneSymbol
from library.enums import ModificationOperation
from pathtests.models import (
    Case,
    PathologyTest,
    PathologyTestGeneModificationOutcome,
    PathologyTestGeneModificationRequest,
    PathologyTestVersion,
    cases_for_user,
)
from patients.models import FollowLeadScientist, Patient
from snpdb.models import ImportStatus


def _gene_list(user, name, *symbols) -> GeneList:
    gene_list = GeneList.objects.create(name=name, user=user, import_status=ImportStatus.SUCCESS)
    for symbol in symbols:
        gene_symbol, _ = GeneSymbol.objects.get_or_create(symbol=symbol)
        GeneListGeneSymbol.objects.create(gene_list=gene_list, gene_symbol=gene_symbol, original_name=symbol)
    return gene_list


def _confirmed_test(curator, name, *symbols) -> PathologyTestVersion:
    """ A test whose v1 is confirmed and active (so its gene list is locked) """
    pathology_test = PathologyTest.objects.create(name=name, curator=curator)
    ptv = PathologyTestVersion.objects.create(pathology_test=pathology_test, gene_list=_gene_list(curator, name, *symbols),
                                              confirmed_date=timezone.now())
    ptv.set_as_active_test()
    return ptv


def _request(ptv, user, symbol, operation) -> PathologyTestGeneModificationRequest:
    gene_symbol, _ = GeneSymbol.objects.get_or_create(symbol=symbol)
    return PathologyTestGeneModificationRequest.objects.create(pathology_test_version=ptv, operation=operation,
                                                               gene_symbol=gene_symbol, user=user)


@override_settings(PATHOLOGY_TESTS_ENABLED=True)
class PathologyTestVersionTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.curator = User.objects.create_user("pathtests_curator")
        cls.requester = User.objects.create_user("pathtests_requester")

    def setUp(self):
        self.client.force_login(self.curator)

    def test_modification_request_str(self):
        ptv = _confirmed_test(self.curator, "str_test", "GATA2")
        gmr = _request(ptv, self.requester, "RUNX1", ModificationOperation.ADD)
        self.assertEqual("str_test (v1) Add RUNX1: Pending", str(gmr))

    def test_test_created_from_confirmed_version_accepts_requests(self):
        source = _confirmed_test(self.curator, "source_test", "GATA2")
        self.assertTrue(source.gene_list.locked)

        response = self.client.post(reverse("manage_pathology_tests"),
                                    {"name": "cloned_test", "pathology_test_version": source.pk})
        self.assertEqual(response.status_code, 302)
        ptv = PathologyTestVersion.objects.get(pathology_test__name="cloned_test")
        self.assertFalse(ptv.gene_list.locked)

        _request(ptv, self.requester, "RUNX1", ModificationOperation.ADD)
        response = self.client.post(ptv.get_absolute_url(), {"add-RUNX1": "accept"})
        self.assertEqual(response.status_code, 200)
        self.assertEqual({"GATA2", "RUNX1"}, set(ptv.gene_list.get_gene_names()))

    def test_delete_then_restore_keeps_latest_confirmed_version_active(self):
        ptv = _confirmed_test(self.curator, "restore_test", "GATA2")
        pathology_test = ptv.pathology_test
        latest_url = reverse("api_view_latest_pathology_test_version", kwargs={"name": pathology_test.name})

        self.client.post(pathology_test.get_absolute_url(), {"delete_test": "true", "delete_text": "delete"})
        pathology_test.refresh_from_db()
        self.assertTrue(pathology_test.deleted)
        self.assertEqual(404, self.client.get(latest_url).status_code)

        self.client.post(pathology_test.get_absolute_url(), {"restore_test": "true"})
        pathology_test.refresh_from_db()
        self.assertFalse(pathology_test.deleted)
        self.assertEqual(ptv, pathology_test.get_active_test_version())
        self.assertEqual(200, self.client.get(latest_url).status_code)

    def test_requests_are_applied_by_operation(self):
        """ A REMOVE for a gene not in the list is a removal (a no-op), never an addition """
        ptv = _confirmed_test(self.curator, "operation_test", "GATA2")
        ptv.confirmed_date = None  # A draft, so requests apply to this version
        ptv.save()
        ptv.gene_list.locked = False
        ptv.gene_list.save()
        _request(ptv, self.requester, "RUNX1", ModificationOperation.REMOVE)

        response = self.client.get(ptv.get_absolute_url())
        self.assertIn("RUNX1", response.context["gene_deletion_requests"])
        self.assertNotIn("RUNX1", response.context["gene_addition_requests"])

        self.client.post(ptv.get_absolute_url(), {"del-RUNX1": "accept"})
        self.assertEqual({"GATA2"}, set(ptv.gene_list.get_gene_names()))

    def test_request_filed_after_page_load_stays_pending(self):
        ptv = _confirmed_test(self.curator, "late_request_test", "GATA2")
        ptv.confirmed_date = None
        ptv.save()
        ptv.gene_list.locked = False
        ptv.gene_list.save()
        answered = _request(ptv, self.requester, "RUNX1", ModificationOperation.ADD)
        late = _request(ptv, self.requester, "PTEN", ModificationOperation.ADD)

        self.client.post(ptv.get_absolute_url(), {"add-RUNX1": "reject"})  # No radio for PTEN in the POST
        answered.refresh_from_db()
        late.refresh_from_db()
        self.assertEqual(PathologyTestGeneModificationOutcome.REJECTED, answered.outcome)
        self.assertEqual(PathologyTestGeneModificationOutcome.PENDING, late.outcome)

    def test_non_curator_post_is_403(self):
        ptv = _confirmed_test(self.curator, "curator_test", "GATA2")
        self.client.force_login(self.requester)
        self.assertEqual(403, self.client.post(ptv.get_absolute_url(), {"confirm_test": "true"}).status_code)
        self.assertEqual(403, self.client.post(ptv.pathology_test.get_absolute_url(),
                                               {"restore_test": "true"}).status_code)


class CasesForUserTest(TestCase):
    def test_cases_led_by_user_or_followed_scientist(self):
        scientist = User.objects.create_user("pathtests_lead_scientist")
        follower = User.objects.create_user("pathtests_follower")
        stranger = User.objects.create_user("pathtests_stranger")
        patient = Patient.objects.create(first_name="Cases", last_name="ForUser")
        case = Case.objects.create(name="led_case", patient=patient, lead_scientist=scientist)

        self.assertEqual([case], list(cases_for_user(scientist)))
        self.assertFalse(cases_for_user(follower).exists())
        FollowLeadScientist.objects.create(user=follower, follow=scientist)
        self.assertEqual([case], list(cases_for_user(follower)))
        self.assertFalse(cases_for_user(stranger).exists())

from django.test import SimpleTestCase, override_settings

from variantgrid.deployment_validation.celery_checks import check_beat_schedule_tasks


class BeatScheduleCheckTest(SimpleTestCase):
    def test_current_schedule_is_valid(self):
        self.assertEqual({"valid": True}, check_beat_schedule_tasks())

    def test_task_module_must_be_imported_by_workers(self):
        """ A task whose module is only reached through urls/views is registered by luck, so it fails the check """
        with override_settings(CELERY_IMPORTS=()):
            check = check_beat_schedule_tasks({
                "heartbeat": {"task": "variantgrid.tasks.server_monitoring_tasks.heartbeat"},
            })
        self.assertFalse(check["valid"])
        self.assertIn("variantgrid.tasks.server_monitoring_tasks is not in settings.CELERY_IMPORTS", check["fix"])

    def test_unknown_and_non_task_names(self):
        check = check_beat_schedule_tasks({
            "missing": {"task": "variantgrid.tasks.server_monitoring_tasks.no_such_task"},
            "not_a_task": {"task": "variantgrid.tasks.server_monitoring_tasks.get_disk_messages"},
        })
        self.assertFalse(check["valid"])
        self.assertIn("no_such_task not found", check["fix"])
        self.assertIn("get_disk_messages is not a celery task", check["fix"])

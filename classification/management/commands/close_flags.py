from django.core.management import BaseCommand
from flags.models import FlagComment, Flag, FlagResolution, FlagStatus, FlagType, FlagTypeResolution


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--user_id', type=int, default=0)
        parser.add_argument('--text', type=str)
        parser.add_argument('--close', action='store_true')
        parser.add_argument('--comment', type=str, default=None)

    def handle(self, *args, **options):
        fcs = FlagComment.objects.filter(text=options["text"])
        flag_qs = Flag.objects.filter(pk__in=fcs.values_list("flag_id", flat=True))
        if user_id := options["user_id"]:
            flag_qs = flag_qs.filter(user_id=user_id)

        open_flags_qs = flag_qs.filter(resolution__status=FlagStatus.OPEN)
        closed_flags_qs = flag_qs.exclude(resolution__status=FlagStatus.OPEN)

        print(f"Open flags found count = {open_flags_qs.count()}")
        print(f"Closed flags found count = {closed_flags_qs.count()}")

        if options["close"]:

            flag_type_to_resolution = {}

            def resolution_for_flag_type(flag_type: FlagType):
                nonlocal flag_type_to_resolution
                if flag_type not in flag_type_to_resolution:
                    resolution = FlagTypeResolution.objects.filter(flag_type=flag_type,
                                                                   resolution__status=FlagStatus.CLOSED).select_related(
                        "resolution").first().resolution
                    flag_type_to_resolution[flag_type] = resolution
                return flag_type_to_resolution[flag_type]

            for flag in open_flags_qs.select_related('flag_type'):
                close_resolution = resolution_for_flag_type(flag.flag_type)
                flag.flag_action(
                    resolution=close_resolution,
                    comment=options["comment"]
                )

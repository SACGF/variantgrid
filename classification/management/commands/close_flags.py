from django.core.management import BaseCommand
from django.db.models import QuerySet

from flags.models import FlagComment, Flag, FlagResolution, FlagStatus, FlagType, FlagTypeResolution


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--user_id', type=int, default=0)
        parser.add_argument('--text', type=str)
        parser.add_argument('--close', action='store_true')
        parser.add_argument('--open', action='store_true')
        parser.add_argument('--comment', type=str, default=None)

    @staticmethod
    def apply_flag_changes(flag_qs: QuerySet[Flag], status: FlagStatus, comment: str):
        flag_type_to_resolution: dict[FlagType, FlagResolution] = {}

        def resolution_for_flag_type(flag_type: FlagType):
            nonlocal flag_type_to_resolution
            if flag_type not in flag_type_to_resolution:
                resolution = FlagTypeResolution.objects.filter(flag_type=flag_type,
                                                               resolution__status=status).select_related(
                    "resolution").first().resolution
                flag_type_to_resolution[flag_type] = resolution
            return flag_type_to_resolution[flag_type]

        for flag in flag_qs.select_related('flag_type'):
            end_resolution = resolution_for_flag_type(flag.flag_type)
            flag.flag_action(
                resolution=end_resolution,
                comment=comment
            )

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
            Command.apply_flag_changes(open_flags_qs, FlagStatus.CLOSED, options["comment"])
        elif options["open"]:
            Command.apply_flag_changes(closed_flags_qs, FlagStatus.OPEN, options["comment"])

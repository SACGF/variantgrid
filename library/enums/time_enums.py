from django.db.models.functions import TruncHour, TruncDay, TruncWeek, TruncMonth, TruncYear


class TimePeriod:
    HOUR = 'h'
    DAY = 'd'
    WEEK = 'w'
    MONTH = 'm'
    YEAR = 'y'

    CHOICES = (
        (HOUR, "Hour"),
        (DAY, "Day"),
        (WEEK, "Week"),
        (MONTH, "Month"),
        (YEAR, "Year"),
    )

    @classmethod
    def truncate_func(cls, time_period):
        TRUNC_OPS = {
            cls.HOUR: TruncHour,
            cls.DAY: TruncDay,
            cls.WEEK: TruncWeek,
            cls.MONTH: TruncMonth,
            cls.YEAR: TruncYear,
        }
        return TRUNC_OPS[time_period]

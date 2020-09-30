from django.contrib.auth.models import User
from library.utils import string_deterministic_hash


class SpaceThemedAvatarProvider:

    @classmethod
    def get_avatar_url(cls, user: User, size):
        if user.username == 'admin_bot':
            return '/static/icons/users/bot.svg'
        icons = [
            "001-alien-of-outer-space.svg",
            "002-astronaut-suit.svg",
            "003-satellite-outer-space-tool-shape-variant.svg",
            "004-small-rocket-ship-silhouette.svg",
            "005-saturn.svg",
            "006-alien.svg",
            "007-ufo.svg",
            "008-alien-1.svg",
            "009-ufo-1.svg",
            "010-alien-2.svg",
            "011-sun.svg",
            "012-rocket.svg",
            "013-spaceship.svg",
            "014-monster.svg",
            "015-monster-1.svg",
            "016-alien-3.svg",
            "017-ufo-2.svg",
            "018-ufo-3.svg",
            "019-alien-4.svg",
            "020-monster-2.svg",
            "021-monster-3.svg",
            "022-alien-5.svg",
            "023-alien-6.svg",
            "024-alien-7.svg",
            "025-alien-8.svg",
            "026-alien-9.svg",
            "027-ufo-4.svg",
            "028-ufo-5.svg",
            "029-alien-10.svg"
        ]
        return '/static/icons/users/' + icons[string_deterministic_hash(user.username) % len(icons)]

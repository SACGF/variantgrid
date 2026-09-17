""" collectstatic, then drop django-compressor's cache entries - one process for the deploy.

    The compressor caches the rendered <link>/<script> tag for each {% compress %} block for
    COMPRESS_REBUILD_TIMEOUT (30 days), naming a bundle under STATIC_ROOT/CACHE, under a key hashed from
    the source files' mtimes. `collectstatic --clear` deletes that directory. With DEBUG off the compressor
    reads the collected files, whose mtimes --clear just refreshed, so the key changes and the entries are
    merely orphaned. With DEBUG on it reads the finders' source files instead, whose mtimes didn't move,
    so every page keeps linking a bundle that 404s until the entries go. Either way they're dead weight
    after collecting, so this drops them.
"""
from compressor.cache import cache as compressor_cache
from compressor.cache import get_cachekey
from django.core.management import CommandParser, call_command
from django.core.management.base import BaseCommand


def clear_compressor_cache() -> int:
    """ Deletes the compressor's keys and nothing else, returning how many went """
    client = compressor_cache._cache.get_client(None, write=True)
    keys = list(client.scan_iter(match=compressor_cache.make_key(get_cachekey("*"))))
    if keys:
        client.delete(*keys)
    return len(keys)


class Command(BaseCommand):
    category = "ops"
    help = "collectstatic, then clear the compressor cache"

    def add_arguments(self, parser: CommandParser):
        parser.add_argument("--noinput", "--no-input", action="store_false", dest="interactive",
                            help="Do NOT prompt the user for input of any kind.")
        parser.add_argument("--clear", action="store_true",
                            help="Clear the existing files using the storage before trying to copy or link the original file.")

    def handle(self, *args, **options):
        call_command("collectstatic", interactive=options["interactive"], clear=options["clear"],
                     verbosity=options["verbosity"])
        print(f"Cleared {clear_compressor_cache()} compressor cache entries")

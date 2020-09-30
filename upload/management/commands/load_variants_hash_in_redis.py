#!/usr/bin/env python3

from django.core.management.base import BaseCommand
import logging

from library.django_utils import get_redis
from upload.tasks.load_variants_hash_in_redis_task import load_variants_hash_in_redis


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--clear', action='store_true')
        parser.add_argument('--async', action='store_true')

    def handle(self, *args, **options):
        root_logger = logging.getLogger('')
        root_logger.setLevel(logging.DEBUG)

        if options["clear"]:
            logging.info("Clearing database")
            redis = get_redis()
            redis.flushall()

        if options["async"]:
            load_variants_hash_in_redis.apply_async()  # @UndefinedVariable
        else:
            load_variants_hash_in_redis()

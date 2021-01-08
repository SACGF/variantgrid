from typing import TypedDict

from django.template import Library
from django.urls import reverse

from snpdb.models import UserSettings, GenomeBuild

register = Library()


@register.inclusion_tag("snpdb/tags/genome_build_url_arg.html")
def genome_build_url_arg(genome_build, url_name, **url_kwargs):
    """ Generates links for user to switch page to their active genome builds
        url_name - must take a parameter 'genome_build_name' """

    builds_with_annotation = GenomeBuild.builds_with_annotation()
    other_genome_builds_exist = builds_with_annotation.exclude(pk=genome_build.pk).exists()

    class BuildUrlDict(TypedDict):
        genome_build: GenomeBuild
        active: bool
        css_class: str
        url: str

    build_urls = []
    for gb in builds_with_annotation:
        url = reverse(url_name, kwargs={"genome_build_name": gb.name, **url_kwargs})
        active = gb == genome_build
        css_class = "active" if active else ""
        build_urls.append({"genome_build": gb, "active": active, "url": url, "css_class": css_class})

    return {
        "build_urls": build_urls,
        "other_genome_builds_exist": other_genome_builds_exist,
    }

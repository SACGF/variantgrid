from django.template import Library

from snpdb.models import UserSettings, GenomeBuild

register = Library()


@register.inclusion_tag("snpdb/tags/genome_build_url_arg.html")
def genome_build_url_arg(user, url_name, genome_build):
    """ Generates links for user to switch page to their active genome builds
        url_name - must take a parameter 'genome_build_name' """

    user_settings = UserSettings.get_for_user(user)
    builds_with_annotation = GenomeBuild.builds_with_annotation()
    other_genome_builds_exist = builds_with_annotation.exclude(pk=genome_build.pk).exists()

    return {
        "url_name": url_name,
        "genome_build": genome_build,
        "builds_with_annotation": builds_with_annotation,
        "user_settings": user_settings,
        "other_genome_builds_exist": other_genome_builds_exist,
    }

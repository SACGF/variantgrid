from django.template import Library

from snpdb.models import UserSettings

register = Library()


@register.inclusion_tag("snpdb/tags/genome_build_url_arg.html")
def genome_build_url_arg(user, url_name, genome_build):
    """ Generates links for user to switch page to their active genome builds
        url_name - must take a parameter 'genome_build_name' """
    user_settings = UserSettings.get_for_user(user)
    has_other_genome_builds = user_settings.get_genome_builds().exclude(pk=genome_build.pk).exists()

    return {
        "url_name": url_name,
        "genome_build": genome_build,
        "user_settings": user_settings,
        "has_other_genome_builds": has_other_genome_builds,
    }

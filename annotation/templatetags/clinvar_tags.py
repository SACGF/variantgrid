from django.template import Library

register = Library()


@register.inclusion_tag("annotation/tags/clinvar_stars.html")
def clinvar_stars(stars):
    MAX_STARS = 4
    stars = ([True] * stars) + ([False] * (MAX_STARS - stars))
    return {"stars": stars}

import uuid

from django.template import Library, loader

from snpdb.models import Wiki

register = Library()


@register.simple_tag(takes_context=True)
def wiki_editor(context, wiki, class_name, unique_keyword, unique_value):
    """ unique_keyword/unique_value must be a way to instantiate a single instance of class_name """

    if not wiki:
        # Load a fake object (won't be saved) so we can check permissions
        klass = Wiki.get_subclass_by_name(class_name)
        wiki = klass(**{unique_keyword: unique_value})

    request = context["request"]
    context = {"uuid": uuid.uuid4(),
               'wiki': wiki,
               'class_name': class_name,
               'unique_keyword': unique_keyword,
               'unique_value': unique_value,
               'has_write_permission': wiki.can_write(request.user)}
    t = loader.get_template(f"snpdb/tags/wiki_tag.html")
    return t.render(context, request=request)

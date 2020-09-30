from django import template


def get_passed_object(token):
    contents = token.split_contents()
    tag_name = contents[0]
    try:
        obj = contents[1]
    except ValueError:
        raise template.TemplateSyntaxError(f"{tag_name} tag requires exactly one argument")

    return obj


def get_passed_objects(token):
    """ Returns a tuple (which you can pass as args) of passed objects """
    contents = token.split_contents()
    return tuple(contents[1:])

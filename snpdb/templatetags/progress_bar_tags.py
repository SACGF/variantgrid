from django.template.library import Library

register = Library()

@register.inclusion_tag("snpdb/tags/progress_bar.html")
def progress_bar(current_step: str, steps: list):
    step_dicts = []
    for step in steps:
        css_classes = ["step"]
        if step == current_step:
            css_classes.append("current")

        step_dicts.append({"css_class": " ".join(css_classes),
                           "name": step})

    return {"steps": step_dicts}

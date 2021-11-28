import markdown
from django import template

from library import tag_utils
from library.tag_utils import get_passed_objects

register = template.Library()


@register.tag
def render_node_text(_parser, token):
    return NodeTextDescriptionTag(tag_utils.get_passed_object(token))


@register.tag
def render_dataframe(_parser, token):
    args = get_passed_objects(token)
    return PandasDataFrameTableTag(*args)


class NodeTextDescriptionTag(template.Node):

    def __init__(self, node):
        self.node_variable = template.Variable(node)

    def render(self, context):
        try:
            node = self.node_variable.resolve(context)
            output_count = node.get_output_count()
            if node.has_input():
                input_count = node.get_input_count()
                if input_count == 0:
                    input_summary = "No inputs"
                    filter_summary = "(No filtering performed on input of 0 variants)"
                else:
                    input_summary = "<UL>Inputs:"
                    for p in node.get_parent_subclasses():
                        input_summary += f"<LI>{p.get_identifier()} ({p.get_output_count()} variants)</LI>"
                    input_summary += "</UL>"
                    filtered = input_count - output_count
                    filter_summary = " (filtered %.1f%% of %d)" % (100.0 * filtered / input_count, input_count)
            else:
                input_summary = ""
                filter_summary = ""

            try:
                node_description_html = markdown.markdown(node.nodewiki.markdown)
            except:
                node_description_html = ''

            name = f"{node.get_class_name()}: {node.name}"
            context = {"node_id": node.pk,
                       "node_name": name,
                       "node_description_html": node_description_html,
                       "node_method": node.get_method_summary(),
                       "filter_summary": filter_summary,
                       "input_summary": input_summary,
                       "output_variants": output_count}

            text = """
                <div node_id="%(node_id)s">
                    <h4>%(node_name)s</h4>
                    <div>%(node_description_html)s</div>

                    <div>%(node_method)s</div>
                    %(input_summary)s
                    <div>Variants: %(output_variants)d%(filter_summary)s
                </div>
            """
            return text % context
        except template.VariableDoesNotExist:
            return ''


class PandasDataFrameTableTag(template.Node):

    def __init__(self, df, significant_figures=None):
        self.df_variable = template.Variable(df)
        if significant_figures is not None:
            self.significant_figures_variable = template.Variable(significant_figures)
        else:
            self.significant_figures_variable = None

    def render(self, context):
        try:
            df = self.df_variable.resolve(context)
            kwargs = {}
            if self.significant_figures_variable is not None:
                significant_figures = self.significant_figures_variable.resolve(context)
                if significant_figures is not None:
                    format_string = "%%.%df" % significant_figures
                    kwargs['float_format'] = lambda f: format_string % f
            return df.to_html(**kwargs)
        except template.VariableDoesNotExist:
            return ''

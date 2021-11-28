import uuid

from django.template import Library

from library.pandas_utils import get_rows_percent_dataframe

register = Library()


@register.inclusion_tag("annotation/tags/version_diff_column_from_to_graph.html")
def version_diff_column_from_to_graph(title, df):
    x = df.columns.tolist()
    y = df.index.tolist()
    perc_df = get_rows_percent_dataframe(df)
    z = perc_df.as_matrix().tolist()
    labels = df.as_matrix().tolist()

    return {'title': title,
            'uuid': uuid.uuid4(),
            'x': x,
            'y': y,
            'z': z,
            'labels': labels}

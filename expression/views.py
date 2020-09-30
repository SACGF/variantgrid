from django.core.exceptions import PermissionDenied
from django.http.response import HttpResponseRedirect
from django.shortcuts import get_object_or_404, render
from django.urls.base import reverse
from expression.forms import ExpressionFileForm
from expression.graphs.volcano_graph import VolcanoGraph
from expression.models import CuffDiffFile
from library.utils import full_class_name
from snpdb.graphs import graphcache


def expression_files(request):
    return render(request, 'expression/expression_files.html')


def expression_graph(request, expression_file_id):
    graph_class_name = full_class_name(VolcanoGraph)
    cached_graph = graphcache.async_graph(graph_class_name, expression_file_id)
    return HttpResponseRedirect(reverse("cached_generated_file_check", kwargs={"cgf_id": cached_graph.id}))


def view_expression_file(request, expression_file_id):
    cuff_diff_file = get_object_or_404(CuffDiffFile, pk=expression_file_id)
    if not request.user.has_perm('view_cuff_diff_file', cuff_diff_file):
        raise PermissionDenied()

    cuff_diff_file_form = ExpressionFileForm(instance=cuff_diff_file)

    context = {'expression_file': cuff_diff_file,
               'expression_file_form': cuff_diff_file_form}
    return render(request, 'expression/view_expression_file.html', context)

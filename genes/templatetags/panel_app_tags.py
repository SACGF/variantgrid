from django.template import Library

from genes.models import PanelAppServer, GeneSymbol

register = Library()


@register.inclusion_tag("genes/tags/panel_app_gene_evidence.html")
def panel_app_gene_evidence(server: PanelAppServer, gene_symbol: GeneSymbol, summary_func: str = None):
    return {
        "server": server,
        "gene_symbol": gene_symbol,
        "summary_func": summary_func,
    }

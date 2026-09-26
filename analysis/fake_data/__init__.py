"""
Fake analyses and variant tags for 'manage.py create_fake_data' (@see library.fake_data): the "analysis" step here,
an analysis over the fake trio, and the "tags" step in variant_tags.py.
"""
from analysis.fake_data import variant_tags  # registers the "tags" step
from analysis.models.enums import SetOperations, TrioInheritance
from analysis.models.models_analysis import Analysis
from analysis.models.nodes.analysis_node import AnalysisNode
from analysis.models.nodes.filters.damage_node import DamageNode
from analysis.models.nodes.filters.population_node import PopulationNode
from analysis.models.nodes.filters.venn_node import VennNode
from analysis.models.nodes.node_utils import update_analysis
from analysis.models.nodes.sources.trio_node import TrioNode
from annotation.models.damage_enums import PathogenicityImpact
from library.fake_data import FakeData, FakeDataContext, register
from snpdb.fake_data import FAKE_USERS

FAKE_ANALYSIS_NAME = "Fake trio analysis"


@register
class FakeAnalysis(FakeData):
    name = "analysis"
    help = "An analysis over the fake trio: dominant and de novo trio nodes, population and impact filters, a venn"
    requires = ("trio",)

    def create(self, context: FakeDataContext, **options):
        trio = context.trio
        if analysis := Analysis.objects.filter(name=FAKE_ANALYSIS_NAME, user=trio.user).first():
            context.stdout.write(f"{analysis} already exists")
            return

        analysis = Analysis(genome_build=context.genome_build, name=FAKE_ANALYSIS_NAME)
        analysis.set_defaults_and_save(trio.user)

        dominant = TrioNode.objects.create(analysis=analysis, trio=trio, inheritance=TrioInheritance.DOMINANT,
                                           x=50, y=50)
        denovo = TrioNode.objects.create(analysis=analysis, trio=trio, inheritance=TrioInheritance.DENOVO,
                                         x=400, y=50)
        rare = _add_child(dominant, PopulationNode(analysis=analysis, percent=1, x=50, y=175))
        damaging = _add_child(rare, DamageNode(analysis=analysis, impact_min=PathogenicityImpact.MODERATE,
                                               x=50, y=300))
        venn = VennNode.objects.create(analysis=analysis, set_operation=SetOperations.UNION, x=225, y=425)
        venn.add_parent(damaging, side=VennNode.LEFT_PARENT)
        venn.add_parent(denovo, side=VennNode.RIGHT_PARENT)
        venn.save()

        update_analysis(analysis.pk)
        context.stdout.write(f"Created {analysis} with {analysis.analysisnode_set.count()} nodes")

    def delete(self, context: FakeDataContext, **options):
        deleted, _ = Analysis.objects.filter(name=FAKE_ANALYSIS_NAME, user__username__in=FAKE_USERS).delete()
        context.stdout.write(f"Deleted the fake analysis ({deleted} rows)")


def _add_child(parent: AnalysisNode, child: AnalysisNode) -> AnalysisNode:
    child.save()
    child.add_parent(parent)
    child.save()
    return child

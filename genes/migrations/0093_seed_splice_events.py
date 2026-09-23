from django.db import migrations

# SACGF/variantgrid#1875 - the splice junctions DRAGEN TSO 500 writes into a CombinedVariantOutput.
# The CombinedVariantOutput keeps passing calls on EGFR, MET and AR only (Illumina's stated rule -
# @see upload.tasks.import_dragen_tso500_combined_variant_output_task), so these are the junctions
# it names; the caller itself runs across the RNA panel. A junction with no row here still imports,
# labelled with its own coordinates (@see genes.gene_splice).
#
# The caller writes breakpoint 1 as the last base of the 5' exon and breakpoint 2 as the 3' exon's
# start exactly as cdot stores it (0-based), so both builds' numbers are exon boundaries of the
# gene's MANE transcript - checked against its exon table, except AR's acceptor, which is the
# cryptic exon CE3 and so is not in one. The intron it sits in is colinear between the builds
# (AR exons 3 and 4 shift by the same 780,158), so the acceptor shifts with them.
SPLICE_EVENTS = [
    # AR-V7: exon 3 spliced to cryptic exon 3 (NM_000044.6 exon 3 ends X:66905968 in GRCh37)
    {"gene_symbol": "AR", "label": "V7", "display": "AR-V7 splice variant", "contig": "X",
     "GRCh37": (66905968, 66914514), "GRCh38": (67686126, 67694672)},
    # EGFRvIII: exon 1 spliced to exon 8, ie exons 2-7 skipped (NM_005228.5)
    {"gene_symbol": "EGFR", "label": "vIII", "display": "EGFRvIII splice variant", "contig": "7",
     "GRCh37": (55087058, 55223522), "GRCh38": (55019365, 55155829)},
    # MET exon 14 skipping: exon 13 spliced to exon 15 (NM_000245.4)
    {"gene_symbol": "MET", "label": "ex14skip", "display": "MET exon 14 skipping", "contig": "7",
     "GRCh37": (116411708, 116414934), "GRCh38": (116771654, 116774880)},
]


def _seed_splice_events(apps, _schema_editor):
    Contig = apps.get_model("snpdb", "Contig")
    GeneSymbol = apps.get_model("genes", "GeneSymbol")
    GenomeBuild = apps.get_model("snpdb", "GenomeBuild")
    SpliceEvent = apps.get_model("genes", "SpliceEvent")

    for data in SPLICE_EVENTS:
        gene_symbol, _ = GeneSymbol.objects.get_or_create(symbol=data["gene_symbol"])
        for genome_build in GenomeBuild.objects.filter(name__in=("GRCh37", "GRCh38")):
            breakpoints = data.get(genome_build.name)
            contig = Contig.objects.filter(genomebuildcontig__genome_build=genome_build,
                                           name=data["contig"]).first()
            if breakpoints is None or contig is None:
                continue
            donor, acceptor = breakpoints
            SpliceEvent.objects.update_or_create(
                genome_build=genome_build, contig=contig, donor=donor, acceptor=acceptor,
                defaults={"gene_symbol": gene_symbol, "label": data["label"],
                          "display": data["display"]})


def _delete_splice_events(apps, _schema_editor):
    SpliceEvent = apps.get_model("genes", "SpliceEvent")
    SpliceEvent.objects.filter(label__in=[d["label"] for d in SPLICE_EVENTS]).delete()


class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0092_spliceevent'),
        ('snpdb', '0263_mandatory_columns_all_collections'),
    ]

    operations = [
        migrations.RunPython(_seed_splice_events, reverse_code=_delete_splice_events),
    ]

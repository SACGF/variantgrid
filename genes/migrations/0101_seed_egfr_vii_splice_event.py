from django.db import migrations

# sapath#457 - EGFRvII, the second EGFR junction Molecular Oncology report by name. Same shape as
# 0093 (breakpoint 1 the last base of the 5' exon, breakpoint 2 the 3' exon's start, both 0-based as
# cdot stores them), with the canonical label 0095 introduced.
#
# Exon boundaries of NM_005228.5 in each build's own alignment - the transcript MO report, and EGFR's
# MANE Select in GRCh38. Exons 14-15 are 249 bases, so the skip is in frame: aa 521-603 of the mature
# protein (545-627 counting from the initiator Met).
SPLICE_EVENTS = [
    # EGFRvII: exon 13 spliced to exon 16, ie exons 14-15 skipped (NM_005228.5)
    {"gene_symbol": "EGFR", "label": "v_ii", "display": "EGFRvII splice variant", "contig": "7",
     "GRCh37": (55229324, 55238867), "GRCh38": (55161631, 55171174)},
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
    for data in SPLICE_EVENTS:
        SpliceEvent.objects.filter(gene_symbol=data["gene_symbol"], label=data["label"]).delete()


class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0100_one_off_transcript_version_modified_cdot_version'),
    ]

    operations = [
        migrations.RunPython(_seed_splice_events, reverse_code=_delete_splice_events),
    ]

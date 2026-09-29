"""
What the fake classification steps share: allele infos resolved onto the fake variants (so records reach the
classification listing and allele pages without ClinGen or HGVS matching), and deleting fake classifications along
with everything built on them.
"""
from contextlib import contextmanager

from django.contrib.auth.models import Permission
from django.contrib.contenttypes.models import ContentType
from django.db.models import QuerySet
from guardian.models import GroupObjectPermission

from annotation.models.models import VariantAnnotation, VariantAnnotationVersion
from classification.models.classification import Classification, ClassificationModification
from classification.models.classification_grouping import (
    AlleleOriginGrouping,
    ClassificationGrouping,
)
from classification.models.classification_reclassification_models import ReclassificationEvent
from classification.models.classification_variant_info_models import (
    ImportedAlleleInfo,
    ImportedAlleleInfoStatus,
    ResolvedVariantInfo,
)
from classification.models.overlaps_model import OverlapContribution
from library.guardian_utils import all_users_group
from snpdb.fake_data import delete_unused_fake_alleles, fake_alleles_for_variants
from snpdb.models.models_genome import GenomeBuild, GenomeBuildPatchVersion

BATCH_SIZE = 2000


def create_fake_allele_infos(genome_build: GenomeBuild,
                             variants_and_genes: list[tuple[int, str]]) -> list[ImportedAlleleInfo]:
    """ An allele info per (variant_id, gene_symbol), already matched to the variant's allele, with the c.HGVS
        from the variant's annotation. The gene chart reads gene_symbol off it, and a record only reaches the
        classification listing once it has an allele to be grouped under """
    variant_ids = [variant_id for variant_id, _ in variants_and_genes]
    allele_ids = fake_alleles_for_variants(genome_build, variant_ids)
    hgvs_c_by_variant = dict(VariantAnnotation.objects.filter(version=VariantAnnotationVersion.latest(genome_build),
                                                              variant_id__in=variant_ids, hgvs_c__isnull=False)
                             .values_list("variant_id", "hgvs_c"))
    patch_version = GenomeBuildPatchVersion.get_or_create(genome_build.name)

    allele_infos = ImportedAlleleInfo.objects.bulk_create([
        ImportedAlleleInfo(imported_c_hgvs=hgvs_c_by_variant.get(variant_id, f"{gene_symbol}:c.{variant_id}A>G"),
                           imported_genome_build_patch_version=patch_version,
                           allele_id=allele_ids[variant_id],
                           status=ImportedAlleleInfoStatus.MATCHED_ALL_BUILDS)
        for variant_id, gene_symbol in variants_and_genes], batch_size=BATCH_SIZE)

    variant_infos = ResolvedVariantInfo.objects.bulk_create([
        ResolvedVariantInfo(allele_info=allele_info, genome_build=genome_build, variant_id=variant_id,
                            gene_symbol_id=gene_symbol, resolved_hgvs=allele_info.imported_c_hgvs)
        for (variant_id, gene_symbol), allele_info in zip(variants_and_genes, allele_infos)], batch_size=BATCH_SIZE)

    build_field = "grch38" if genome_build.name == "GRCh38" else "grch37"
    for allele_info, variant_info in zip(allele_infos, variant_infos):
        setattr(allele_info, build_field, variant_info)
    ImportedAlleleInfo.objects.bulk_update(allele_infos, [build_field], batch_size=BATCH_SIZE)
    return allele_infos


def assign_all_users_permissions(klass, object_ids: list[int]) -> int:
    """ bulk_create skips the permissions mixin's save(), and everyone can see (and clean up) fake data """
    group = all_users_group()
    content_type = ContentType.objects.get_for_model(klass)
    permissions = Permission.objects.filter(content_type=content_type,
                                            codename__in=[klass.get_read_perm(), klass.get_write_perm()])
    total = 0
    for permission in permissions:
        for i in range(0, len(object_ids), BATCH_SIZE):
            batch = object_ids[i:i + BATCH_SIZE]
            GroupObjectPermission.objects.bulk_create(
                [GroupObjectPermission(group=group, permission=permission, content_type=content_type,
                                       object_pk=str(pk)) for pk in batch])
            total += len(batch)
    return total


def delete_fake_classifications(classifications_qs: QuerySet[Classification]) -> str:
    """ Deletes the classifications and what was built on them - permissions, timelines, allele infos, empty
        groupings, and the alleles once no other fake data uses them. Returns a summary """
    classification_ids = list(classifications_qs.values_list("pk", flat=True))
    modification_ids = list(ClassificationModification.objects
                            .filter(classification__in=classification_ids).values_list("pk", flat=True))
    allele_info_ids = [pk for pk in classifications_qs.values_list("allele_info_id", flat=True) if pk]
    allele_ids = [pk for pk in ImportedAlleleInfo.objects.filter(pk__in=allele_info_ids)
                  .values_list("allele_id", flat=True) if pk]
    lab_ids = set(classifications_qs.values_list("lab_id", flat=True))

    deleted_permissions = 0
    for klass, object_ids in ((Classification, classification_ids),
                              (ClassificationModification, modification_ids)):
        content_type = ContentType.objects.get_for_model(klass)
        object_pks = [str(pk) for pk in object_ids]
        for i in range(0, len(object_pks), BATCH_SIZE):
            deleted, _ = GroupObjectPermission.objects.filter(
                content_type=content_type, object_pk__in=object_pks[i:i + BATCH_SIZE]).delete()
            deleted_permissions += deleted

    ReclassificationEvent.objects.filter(classification__in=classification_ids).delete()
    Classification.objects.filter(pk__in=classification_ids).delete()  # cascades to modifications
    ResolvedVariantInfo.objects.filter(allele_info_id__in=allele_info_ids).delete()
    ImportedAlleleInfo.objects.filter(pk__in=allele_info_ids).delete()
    empty_groupings_qs = ClassificationGrouping.objects.filter(
        lab__in=lab_ids, allele_origin_grouping__allele__in=allele_ids, classificationgroupingentry__isnull=True)
    # while the grouping still says which lab they're from - the audit log names them on delete
    OverlapContribution.objects.filter(classification_grouping__in=empty_groupings_qs).delete()
    empty_groupings_qs.delete()
    AlleleOriginGrouping.objects.filter(allele__in=allele_ids, classificationgrouping__isnull=True).delete()
    deleted_alleles = delete_unused_fake_alleles(allele_ids)
    return (f"Deleted {len(classification_ids)} classifications, {len(modification_ids)} modifications, "
            f"{deleted_permissions} permissions, {len(allele_info_ids)} allele infos and {deleted_alleles} alleles")


@contextmanager
def keeping_our_timestamps(klass):
    """ created/modified are auto_now_add/auto_now, which would stamp years of curation as happening now """
    created = klass._meta.get_field("created")
    modified = klass._meta.get_field("modified")
    created.auto_now_add = False
    modified.auto_now = False
    try:
        yield
    finally:
        created.auto_now_add = True
        modified.auto_now = True

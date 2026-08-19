import logging
from collections import Counter, defaultdict
from functools import cached_property
from typing import Optional

from django.conf import settings
from django.db.models import Q, StringAgg, TextField, Value
from django.urls import reverse

from annotation.models.models import VariantAnnotation
from annotation.models.models_phenotype_match import PATIENT_ONTOLOGY_TERM_PATH
from classification.models import ClassificationModification
from classification.models.evidence_mixin_summary_cache import clinical_significance_pills
from library.log_utils import log_traceback
from library.unit_percent import format_af
from ontology.models import OntologyService
from patients.models import Patient
from patients.models_enums import NucleicAcid, TissueStatus, Zygosity
from snpdb.models import (
    CohortGenotype,
    CohortGenotypeCollection,
    CohortSample,
    Locus,
    Sample,
    Variant,
    VariantZygosityCount,
    VariantZygosityCountCollection,
    VCFFilter,
)
from upload.models import ModifiedImportedVariant
from variantgrid.perm_path import get_visible_url_names

SAMPLE_ENRICHMENT_KIT_PATH = "samplefromsequencingsample__sequencing_sample__enrichment_kit__name"
# Sample fields a row copies across under a "sample__" prefix
COPY_SAMPLE_FIELDS = ["id", "name", "patient", "patient__patient_code", SAMPLE_ENRICHMENT_KIT_PATH]
# Patient / specimen / extraction identifiers behind the grid's expandable row detail (private#2837).
# Copied across under their own names, and only shown for patients the user can view
SAMPLE_DETAIL_FIELDS = [
    "patient__external_pk__code",
    "extraction",
    "extraction__reference_id",
    "extraction__external_pk__code",
    "extraction__nucleic_acid_source",
    "extraction__specimen",
    "extraction__specimen__patient",
    "extraction__specimen__patient__patient_code",
    "extraction__specimen__patient__external_pk__code",
    "extraction__specimen__reference_id",
    "extraction__specimen__external_pk__code",
    "extraction__specimen__tissue__name",
    "extraction__specimen__tissue_status",
    "extraction__specimen__collection_date",
]
# Where a row's patient comes from, in order of preference - the sequencing pipeline links a sample to
# an Extraction without touching Sample.patient, so the specimen's patient is often the only one there
PATIENT_DETAIL_PATHS = [
    ("sample__patient", "sample__patient__patient_code", "patient__external_pk__code"),
    ("extraction__specimen__patient", "extraction__specimen__patient__patient_code",
     "extraction__specimen__patient__external_pk__code"),
]
# CohortGenotype rows are streamed in chunks, so an early stop leaves the rest of a common variant's
# VCFs unfetched. Collections are bulk loaded per chunk, so this also sets the queries-per-scan rate.
COHORT_GENOTYPE_CHUNK_SIZE = 100
# Zygosity.CHOICES uses '?' for unknown, which doesn't make a great JSON key or column header
LOCUS_COUNT_ZYGOSITY_LABELS = [(Zygosity.HOM_REF, "REF"), (Zygosity.HET, "HET"),
                               (Zygosity.HOM_ALT, "HOM_ALT"), (Zygosity.UNKNOWN_ZYGOSITY, "Unknown")]


class VariantZygosityCounts:
    """ Permission-scoped per-variant sample/zygosity lookup - the reusable core extracted
        from VariantSampleGenotypes (no presentation concerns). Shared by the
        Variantopedia sample-information API (subclass below) and the Beacon
        g_variants endpoint (§5.2 stage 2). Gates by Sample.filter_for_user (anonymous ->
        public group), so presence/counts are only ever asserted from samples the requester
        may read. """

    def __init__(self, user, variant, max_rows: Optional[int] = None, stop_when_full: bool = False):
        """ max_rows caps how many visible rows are *materialised* - None for all, 0 for counts only.
            Counting reads the packed genotype arrays directly, so it never depends on max_rows.

            stop_when_full also stops the CohortGenotype scan once max_rows is reached, which leaves
            the counts partial (self.partial) - the caller has to source its totals elsewhere. That's
            the samples grid, where a common variant means scanning thousands of collections nobody
            is going to look at. Callers that need exact counts (Beacon) leave it False. """
        self.user = user
        self.variant = variant
        # Builds can share contigs (GRCh37/38 share MT and some unplaced scaffolds) in which case both builds
        # have the same variant. Scope to every build the variant is in - scoping to one would count the other
        # build's samples as observations without ever making them visible rows, reporting them as hidden.
        self.genome_builds = sorted(variant.genome_builds, key=lambda gb: gb.name)
        self.num_samples = Sample.objects.filter(vcf__genome_build__in=self.genome_builds).count()
        visible_samples_qs = Sample.filter_for_user(user).filter(vcf__genome_build__in=self.genome_builds)
        self.user_sample_ids = set(visible_samples_qs.values_list("pk", flat=True))
        self.num_user_samples = len(self.user_sample_ids)

        self.visible_zygosity_counter = Counter()  # raw {zygosity: count} for visible samples
        self.locus_counter = defaultdict(Counter)
        self.num_observations = 0
        self.num_visible_observations = 0
        self.partial = False  # Stopped scanning early, so the counts above are of what we scanned
        self._collections = {}  # cgc_id -> CohortGenotypeCollection
        self._packed_samples = defaultdict(list)  # cgc_id -> [(sample_id, packed index)], row order

        rows_to_build = []
        for cg_values_chunk in self._iter_cohort_genotype_chunks(variant):
            self._load_collections({cg_values["variant__cohortgenotype__collection"]
                                    for cg_values in cg_values_chunk})
            for cg_values in cg_values_chunk:
                variant_id = cg_values["variant"]
                cgc = self._collections.get(cg_values["variant__cohortgenotype__collection"])
                if cgc is not None and self._cohort_out_of_date(cgc):
                    continue

                # Alleles at the locus with no CohortGenotype have no samples - the defaultdict touch is what
                # gives them a locus counts row, with zeros rather than being counted as unknown zygosity
                zygosity_counter = self.locus_counter[variant_id]
                if cgc is None:
                    continue

                samples_zygosity = cg_values["variant__cohortgenotype__samples_zygosity"]
                for sample_id, index in self._packed_samples[cgc.pk]:
                    zygosity = samples_zygosity[index]
                    if zygosity == Zygosity.MISSING:
                        continue
                    zygosity_counter[zygosity] += 1
                    if variant_id != variant.pk:
                        continue

                    self.num_observations += 1
                    if sample_id in self.user_sample_ids:
                        self.visible_zygosity_counter[zygosity] += 1
                        self.num_visible_observations += 1
                        if max_rows is None or len(rows_to_build) < max_rows:
                            rows_to_build.append((cgc, cg_values, sample_id, index, zygosity))
                        elif stop_when_full:
                            self.partial = True
                            break
                if self.partial:
                    break
            if self.partial:
                break

        self.visible_rows = self._build_rows(rows_to_build)

    @cached_property
    def genome_builds_str(self) -> str:
        return ", ".join(str(gb) for gb in self.genome_builds)

    @property
    def has_hidden_samples(self) -> bool:
        """ Are there samples in these builds the user can't see (regardless of this variant) """
        return self.num_samples > self.num_user_samples

    @property
    def truncated(self) -> bool:
        """ Were there visible observations we didn't materialise a row for (max_rows)? """
        return self.partial or len(self.visible_rows) < self.num_visible_observations

    @cached_property
    def global_locus_counts(self) -> dict[int, VariantZygosityCount]:
        """ Precomputed ref/het/hom/unk per variant at this locus, keyed on variant id.

            These are global - every sample in the database, no permission scoping - which is what the
            locus counts have always been (the scan counts every sample at the locus, and only filters
            rows). So when we stop scanning early these are the same numbers, for one indexed query. """
        try:
            vzcc = VariantZygosityCountCollection.get_global_germline_counts()
        except VariantZygosityCountCollection.DoesNotExist:
            return {}
        vzc_qs = VariantZygosityCount.objects.filter(collection=vzcc, variant__locus=self.variant.locus)
        return {vzc.variant_id: vzc for vzc in vzc_qs}

    @property
    def global_num_observations(self) -> Optional[int]:
        """ Every call for this variant in the database, or None if the counts aren't available """
        if vzc := self.global_locus_counts.get(self.variant.pk):
            return vzc.ref_count + vzc.het_count + vzc.hom_count + vzc.unk_count
        return None

    @property
    def num_invisible_observations(self) -> int:
        return self.num_observations - self.num_visible_observations

    @property
    def num_visible_het(self) -> int:
        return self.visible_zygosity_counter[Zygosity.HET]

    @property
    def num_visible_hom_alt(self) -> int:
        return self.visible_zygosity_counter[Zygosity.HOM_ALT]

    @property
    def num_visible_alt(self) -> int:
        """ Visible het + hom_alt observations - the presence/count a Beacon reports. """
        return self.num_visible_het + self.num_visible_hom_alt

    @property
    def has_visible_observations(self) -> bool:
        return self.num_visible_alt > 0

    @property
    def has_observations(self) -> bool:
        return self.num_observations > 0

    @staticmethod
    def _get_sample_values_for_variant_via_cohort_genotype(locus_qs):
        """ This is the new, preferred way - as it gets all of the samples for a VCF at once.

            Archived VCFs naturally fall out: archive deletes the CohortGenotypeCollection rows
            (and partitions), so the cohortgenotype join yields no rows for archived VCFs. """
        no_cohort = Q(variant__cohortgenotype__isnull=True)
        vcf_cohort = Q(variant__cohortgenotype__collection__cohort__vcf__isnull=False)
        qs = locus_qs.filter(no_cohort | vcf_cohort)  # Only VCF CohortGenotypes (not generated cohorts)
        return qs.values("variant",
                         "variant__alt",
                         "variant__cohortgenotype__collection",
                         "variant__cohortgenotype__collection__cohort__vcf__name",
                         "variant__cohortgenotype__collection__cohort__vcf__genome_build",
                         "variant__cohortgenotype__filters",
                         "variant__cohortgenotype__samples_filters",
                         "variant__cohortgenotype__samples_zygosity",
                         "variant__cohortgenotype__samples_allele_frequency",
                         "variant__cohortgenotype__samples_allele_depth",
                         "variant__cohortgenotype__samples_read_depth",
                         "variant__cohortgenotype__samples_phred_likelihood")

    def _iter_cohort_genotype_chunks(self, variant):
        """ Streamed, so stopping early doesn't fetch the packed arrays of every remaining VCF """
        locus_qs = Locus.objects.filter(pk=variant.locus.pk)
        values_qs = self._get_sample_values_for_variant_via_cohort_genotype(locus_qs)
        chunk = []
        for cg_values in values_qs.iterator(chunk_size=COHORT_GENOTYPE_CHUNK_SIZE):
            chunk.append(cg_values)
            if len(chunk) == COHORT_GENOTYPE_CHUNK_SIZE:
                yield chunk
                chunk = []
        if chunk:
            yield chunk

    def _load_collections(self, cgc_ids: set):
        """ Everything counting needs, in a fixed number of queries however many VCFs are in the chunk:
            the collections, and per collection the (sample_id, packed genotype index) pairs in sample
            pk order. No DNA controls are dropped here, as they were never counted. """
        # None is a variant at the locus nobody has; the rest can repeat across chunks
        cgc_ids = {cgc_id for cgc_id in cgc_ids if cgc_id is not None and cgc_id not in self._collections}
        if not cgc_ids:
            return

        cgc_ids_by_cohort_id = defaultdict(list)  # A cohort can have more than one collection (versions)
        for cgc in CohortGenotypeCollection.objects.filter(pk__in=cgc_ids).select_related("cohort", "cohort__vcf"):
            self._collections[cgc.pk] = cgc
            cgc_ids_by_cohort_id[cgc.cohort_id].append(cgc.pk)

        cohort_samples_qs = CohortSample.objects.filter(cohort_id__in=cgc_ids_by_cohort_id,
                                                        sample__no_dna_control=False).order_by("sample_id")
        for cohort_id, sample_id, index in cohort_samples_qs.values_list("cohort_id", "sample_id",
                                                                        "cohort_genotype_packed_field_index"):
            for cgc_id in cgc_ids_by_cohort_id[cohort_id]:
                self._packed_samples[cgc_id].append((sample_id, index))

    @staticmethod
    def _cohort_out_of_date(cgc: CohortGenotypeCollection) -> bool:
        if cgc.cohort_version == cgc.cohort.version:
            return False
        logging.warning("CohortGenotypeCollection (without VCF) %s v.%d out of date with Cohort v.%d",
                        cgc.pk, cgc.cohort_version, cgc.cohort.version)
        # A cohort with a VCF can't change, so that's not really out of date, just a bug - see issue #1053
        return cgc.cohort.vcf is None

    def _build_rows(self, rows_to_build: list[tuple]) -> list[dict]:
        """ The expensive half - patient ontology terms and formatted values, for the rows we'll show """
        sample_values_by_id = self._get_sample_values({sample_id for _cgc, _cg_values, sample_id, _index, _zygosity
                                                       in rows_to_build})
        # Filter codes are per-VCF, so they can't share a lookup - but they can share a query
        filter_code_lookup_by_vcf = defaultdict(dict)
        vcf_ids = {cgc.cohort.vcf_id for cgc, _cg_values, _sample_id, _index, _zygosity in rows_to_build}
        for vcf_id, filter_code, filter_id in VCFFilter.objects.filter(vcf_id__in=vcf_ids) \
                .values_list("vcf_id", "filter_code", "filter_id"):
            filter_code_lookup_by_vcf[vcf_id][filter_code] = filter_id

        rows = []
        for cgc, cg_values, sample_id, index, zygosity in rows_to_build:
            rows.append(self._build_row(cg_values, sample_values_by_id[sample_id], index, zygosity,
                                        filter_code_lookup_by_vcf[cgc.cohort.vcf_id]))
        return rows

    @staticmethod
    def _get_sample_values(sample_ids: set[int]) -> dict[int, dict]:
        """ One query for every row we're building, whichever VCFs they came from """
        if not sample_ids:
            return {}

        ontology_path = f"patient__{PATIENT_ONTOLOGY_TERM_PATH}__name"
        q_hpo = Q(**{f"patient__{PATIENT_ONTOLOGY_TERM_PATH}__ontology_service": OntologyService.HPO})
        q_omim = Q(**{f"patient__{PATIENT_ONTOLOGY_TERM_PATH}__ontology_service": OntologyService.OMIM})
        q_mondo = Q(**{f"patient__{PATIENT_ONTOLOGY_TERM_PATH}__ontology_service": OntologyService.MONDO})
        annotation_kwargs = {"patient_hpo": StringAgg(ontology_path, Value('|'), filter=q_hpo,
                                                      distinct=True, output_field=TextField()),
                             "patient_omim": StringAgg(ontology_path, Value('|'), filter=q_omim,
                                                       distinct=True, output_field=TextField()),
                             "patient_mondo": StringAgg(ontology_path, Value('|'), filter=q_mondo,
                                                        distinct=True, output_field=TextField())}
        samples_qs = Sample.objects.filter(pk__in=sample_ids).order_by("pk").annotate(**annotation_kwargs)
        sample_values = samples_qs.values("vcf__allele_frequency_percent", *COPY_SAMPLE_FIELDS,
                                          *SAMPLE_DETAIL_FIELDS, *list(annotation_kwargs.keys()))
        return {s_values["id"]: s_values for s_values in sample_values}

    @staticmethod
    def _packed_value(cg_values: dict, field: str, index: int):
        if values := cg_values["variant__cohortgenotype__" + field]:
            return values[index]
        return -1

    @classmethod
    def _build_row(cls, cg_values: dict, s_values: dict, index: int, zygosity: str, filter_code_lookup) -> dict:
        # The displayed AF may be percent formatted text, so keep the unit value for graphing
        packed_af = cls._packed_value(cg_values, "samples_allele_frequency", index)
        source_in_percent = s_values["vcf__allele_frequency_percent"]
        allele_frequency_unit = format_af(packed_af, source_in_percent=source_in_percent,
                                          dest_in_percent=False,
                                          missing_value=CohortGenotype.MISSING_NUMBER_VALUE)
        allele_frequency = format_af(packed_af, source_in_percent=source_in_percent,
                                     dest_in_percent=settings.VARIANT_ALLELE_FREQUENCY_CLIENT_SIDE_PERCENT,
                                     missing_value=CohortGenotype.MISSING_NUMBER_VALUE)

        sample_filters = None  # VCF had no FT field
        if samples_filters := cg_values["variant__cohortgenotype__samples_filters"]:
            if (ft := samples_filters[index]) != CohortGenotype.MISSING_FT_VALUE:
                sample_filters = VCFFilter.format_filter_codes(filter_code_lookup, ft)

        row = {
            "variant": cg_values["variant"],
            "genome_build": cg_values["variant__cohortgenotype__collection__cohort__vcf__genome_build"],
            "zygosity": zygosity,
            "allele_frequency": allele_frequency,
            "allele_frequency_unit": allele_frequency_unit,
            "allele_depth": cls._packed_value(cg_values, "samples_allele_depth", index),
            "read_depth": cls._packed_value(cg_values, "samples_read_depth", index),
            "phred_likelihood": cls._packed_value(cg_values, "samples_phred_likelihood", index),
            "filters": VCFFilter.format_filter_codes(filter_code_lookup,
                                                     cg_values["variant__cohortgenotype__filters"]),
            "sample_filters": sample_filters,
            "sample": s_values["id"],
            "sample__vcf__name": cg_values["variant__cohortgenotype__collection__cohort__vcf__name"],
        }
        for k, v in s_values.items():
            if k in COPY_SAMPLE_FIELDS:
                k = "sample__" + k
            row[k] = v
        return row


class VariantSampleGenotypes(VariantZygosityCounts):
    """ Builds the JSON payload for the variant details samples grid - one response per variant.
        The grid is client side, so everything it draws (rows, locus counts, multi-allelic
        variants, sample/observation summary) comes from here. """

    def to_json(self) -> dict:
        rows_json = [self._row_to_json(row) for row in self.visible_rows]
        self._add_classifications(rows_json)

        return {
            "genome_builds": [gb.name for gb in self.genome_builds],
            "variant_id": self.variant.pk,
            "variant": str(self.variant),
            "variant_url": reverse("view_variant", kwargs={"variant_id": self.variant.pk}),
            "g_hgvs": self._get_g_hgvs(),
            "samples": {
                "total": self.num_samples,
                "visible": self.num_user_samples,
            },
            "observations": self._get_observations(),
            "zygosity_counts": self._get_zygosity_counts(),
            "locus_counts": self._get_locus_counts(),
            "multiallelic": self._get_multiallelic(),
            "partial": self.partial,
            "truncated": self.truncated,
            "rows": rows_json,
        }

    def _get_zygosity_counts(self) -> dict:
        """ These label the grid's zygosity filter, so after a partial scan they count the rows we
            actually have rather than the one-past-the-cap observation that stopped us """
        if self.partial:
            return dict(Counter(row["zygosity"] for row in self.visible_rows))
        return dict(self.visible_zygosity_counter)

    def _get_observations(self) -> dict:
        """ A partial scan counted only what it scanned, so the totals come from the precomputed global
            counts instead - and we say nothing about how many of them this user can see, because
            working that out is the scan we just skipped. """
        if not self.partial:
            return {
                "total": self.num_observations,
                "visible": self.num_visible_observations,
                "invisible": self.num_invisible_observations,
            }
        return {
            "total": self.global_num_observations,
            "visible": None,
            "invisible": None,
        }

    def _add_classifications(self, rows_json: list[dict]):
        """ Only ever sees rows we materialised, so the cap bounds this lookup too """
        classifications_by_sample_id = self._get_classifications_by_sample_id({row["sample"] for row in rows_json})
        for row in rows_json:
            row["classifications"] = classifications_by_sample_id.get(row["sample"], [])

    def _get_classifications_by_sample_id(self, sample_ids: set[int]) -> dict[int, list[dict]]:
        """ Scoped to the allele, so classifications curated against another build's variant come through too.
            Classifications not linked to a sample are shown by the classifications section of the page. """
        allele = self.variant.allele
        if allele is None or not sample_ids:
            return {}

        classifications_by_sample_id = defaultdict(list)
        qs = ClassificationModification.latest_for_user(self.user, allele=allele, published=True,
                                                        classification__sample__in=sample_ids)
        for cm in qs.select_related("classification__lab"):
            classification = cm.classification
            pills = clinical_significance_pills(classification.summary_typed, classification.allele_origin_bucket)
            classification_json = {
                "id": classification.pk,
                "url": classification.get_absolute_url(),
                "lab": str(classification.lab),
            }
            classifications_by_sample_id[classification.sample_id].extend(classification_json | pill for pill in pills)
        return classifications_by_sample_id

    def _get_g_hgvs(self) -> Optional[str]:
        """ Only used to name the CSV export - reference variants have to generate it, which needs
            annotation the deployment may not have, and that's not worth losing the samples over """
        try:
            return VariantAnnotation.get_hgvs_g(self.variant)
        except Exception:  # pylint: disable=broad-except
            log_traceback()
            return None

    @cached_property
    def _visible_patient_ids(self) -> set[int]:
        """ A sample carries its VCF's permissions while a patient carries its own, so the identifying
            detail is only given for patients the user may view (private#2837) """
        patient_ids = {row["sample__patient"] for row in self.visible_rows}
        patient_ids |= {row["extraction__specimen__patient"] for row in self.visible_rows}
        patient_ids.discard(None)
        if not patient_ids:
            return set()
        visible_qs = Patient.filter_for_user(self.user).filter(pk__in=patient_ids)
        return set(visible_qs.values_list("pk", flat=True))

    @staticmethod
    def _url_if_visible(url_name: str, **kwargs) -> Optional[str]:
        """ A deployment without patients (eg Shariant) unregisters these urls entirely """
        if get_visible_url_names().get(url_name):
            return reverse(url_name, kwargs=kwargs)
        return None

    def _get_patient_detail(self, row: dict) -> dict:
        for patient_path, code_path, external_pk_path in PATIENT_DETAIL_PATHS:
            if (patient_id := row[patient_path]) in self._visible_patient_ids:
                return {
                    "patient_code": row[code_path] or row[external_pk_path],
                    "patient_url": self._url_if_visible("view_patient", patient_id=patient_id),
                }
        return {}

    def _get_sample_details(self, row: dict) -> dict:
        """ Patient/specimen/extraction identifiers, drawn in the grid's expandable row detail """
        details = self._get_patient_detail(row)
        if row["extraction__specimen__patient"] in self._visible_patient_ids:
            details["specimen"] = (row["extraction__specimen__reference_id"]
                                   or row["extraction__specimen__external_pk__code"])
            details["specimen_url"] = self._url_if_visible("view_specimen",
                                                           specimen_id=row["extraction__specimen"])
            details["specimen_tissue"] = row["extraction__specimen__tissue__name"]
            # Unknown is the default for a specimen nobody has said anything about, so it's not worth a line
            tissue_status = row["extraction__specimen__tissue_status"]
            if tissue_status and tissue_status != TissueStatus.UNKNOWN:
                details["specimen_tissue_status"] = TissueStatus(tissue_status).label
            if collection_date := row["extraction__specimen__collection_date"]:
                details["specimen_collection_date"] = collection_date.strftime(settings.DATE_FORMAT)
            details["extraction"] = row["extraction__reference_id"] or row["extraction__external_pk__code"]
            details["extraction_url"] = self._url_if_visible("view_extraction", extraction_id=row["extraction"])
            if nucleic_acid := row["extraction__nucleic_acid_source"]:
                details["nucleic_acid"] = NucleicAcid(nucleic_acid).label
        return {k: v for k, v in details.items() if v}

    def _row_to_json(self, row: dict) -> dict:
        sample_id = row["sample"]
        allele_frequency_unit = row["allele_frequency_unit"]
        if allele_frequency_unit == CohortGenotype.MISSING_NUMBER_VALUE:
            allele_frequency_unit = None
        return {
            "genome_build": row["genome_build"],
            "sample": sample_id,
            "sample_name": row["sample__name"],
            "sample_url": reverse("view_sample", kwargs={"sample_id": sample_id}),
            "patient": row["sample__patient"],  # Phenotype graphs count distinct patients, not samples
            "vcf": row["sample__vcf__name"],
            "zygosity": row["zygosity"],
            "allele_frequency": row["allele_frequency"],
            "allele_frequency_unit": allele_frequency_unit,  # Graphs need a number, not the display value
            "allele_depth": row["allele_depth"],
            "read_depth": row["read_depth"],
            "phred_likelihood": row["phred_likelihood"],
            "filters": row["filters"],
            "sample_filters": row["sample_filters"],
            "enrichment_kit": row["sample__" + SAMPLE_ENRICHMENT_KIT_PATH],
            "patient_hpo": row["patient_hpo"],
            "patient_omim": row["patient_omim"],
            "patient_mondo": row["patient_mondo"],
        } | self._get_sample_details(row)

    def _get_locus_counts(self) -> list[dict]:
        """ Zygosity counts for every variant at this locus, this variant first """
        counts_by_variant_id = self._get_locus_zygosity_counts()
        # str(v) and v.alt.seq below reach through to the locus/sequence rows
        variant_qs = Variant.objects.filter(pk__in=counts_by_variant_id) \
            .select_related("locus__contig", "locus__ref", "alt")
        variant_by_id = {v.pk: v for v in variant_qs}

        sorted_rows = []
        for variant_id, zygosity_counts in counts_by_variant_id.items():
            v = variant_by_id[variant_id]
            row = {
                "variant_id": variant_id,
                "variant": str(v),
                "url": reverse("view_variant", kwargs={"variant_id": variant_id}),
                "description": "",
                "total": sum(zygosity_counts.values()),
            }
            for zygosity, label in LOCUS_COUNT_ZYGOSITY_LABELS:
                row[label] = zygosity_counts[zygosity]

            if variant_id == self.variant.pk:
                row["description"] = "This variant"
                sort_order = 1
            elif v.is_reference:
                sort_order = 0
            else:
                sort_order = sum(map(ord, v.alt.seq))
            sorted_rows.append((sort_order, row))

        sorted_rows.sort(key=lambda sort_order_row: sort_order_row[0])
        return [row for _, row in sorted_rows]

    def _get_locus_zygosity_counts(self) -> dict[int, Counter]:
        """ From the scan, unless we stopped early - the precomputed global counts are the same
            numbers (the locus counts were never permission scoped), so fall back to those """
        if not self.partial:
            return self.locus_counter

        counts_by_variant_id = {}
        for variant_id, vzc in self.global_locus_counts.items():
            counts_by_variant_id[variant_id] = Counter({Zygosity.HOM_REF: vzc.ref_count,
                                                        Zygosity.HET: vzc.het_count,
                                                        Zygosity.HOM_ALT: vzc.hom_count,
                                                        Zygosity.UNKNOWN_ZYGOSITY: vzc.unk_count})
        return counts_by_variant_id or self.locus_counter

    def _get_multiallelic(self) -> list[dict]:
        """ Variants that were split off this one's VCF record onto another locus """
        by_multiallelic = ModifiedImportedVariant.get_other_loci_variants_by_multiallelic(self.variant)
        multiallelic = []
        for old_multiallelic, variants in by_multiallelic.items():
            variant_list = [{"id": v.pk, "label": str(v), "url": v.get_absolute_url()}
                            for v in sorted(variants, key=lambda v: v.pk)]
            multiallelic.append({"multiallelic": old_multiallelic, "variants": variant_list})
        return multiallelic

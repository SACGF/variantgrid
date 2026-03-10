# Plan: Integrate hgvs_shim into VariantGrid

## Overview

`hgvs_shim` is a library that extracts the pyhgvs/biocommons HGVS wrapper code from VariantGrid
into a standalone, Django-independent package. This plan covers replacing VariantGrid's duplicated
code with imports from `hgvs_shim` while keeping the Django-specific adapter layer in VariantGrid.

hgvs_shim is available on PyPI. Use `hgvs_shim==0.2.0`.

## What hgvs_shim Provides vs What VG Has

### Duplicated code to replace:

| VG Location | hgvs_shim | Action |
|---|---|---|
| `genes/hgvs/hgvs.py`: `HGVSVariant` | `hgvs_shim.HGVSVariant` | Import from hgvs_shim |
| `genes/hgvs/hgvs.py`: `HGVSException*` | `hgvs_shim.HGVSException*` | Import from hgvs_shim |
| `genes/hgvs/pyhgvs/`: `PyHGVSVariant` | `hgvs_shim.PyHGVSVariant` | Import from hgvs_shim |
| `genes/hgvs/biocommons_hgvs/`: `BioCommonsHGVSVariant` | `hgvs_shim.BioCommonsHGVSVariant` | Import from hgvs_shim |
| `genes/hgvs/pyhgvs/`: `PyHGVSConverter` | `hgvs_shim.PyHGVSConverter` | VG wraps via composition |
| `genes/hgvs/biocommons_hgvs/`: `BioCommonsHGVSConverter` | `hgvs_shim.BioCommonsHGVSConverter` | VG subclasses |
| `genes/hgvs/hgvs_converter_combo.py` | `hgvs_shim.ComboCheckerHGVSConverter` | Keep VG version (different API) |

### Key API differences:

1. **`HGVSConverter` method signatures**: hgvs_shim uses `(chrom, position, ref, alt)` and
   `TranscriptInfo` dataclass. VG uses `VariantCoordinate` and `TranscriptVersion` Django models.

2. **Constructor parameters**: hgvs_shim's `PyHGVSConverter(fasta_file, get_transcript)` and
   `BioCommonsHGVSConverter(assembly_name, hdp)`. VG's take `genome_build: GenomeBuild`.

3. **`HGVSConverterType`**: hgvs_shim has PYHGVS=1, BIOCOMMONS_HGVS=2, COMBO=3 as `IntEnum`.
   VG adds `CLINGEN_ALLELE_REGISTRY = 4`.

4. **Return type of HGVS→coordinate**: hgvs_shim returns `(chrom, position, ref, alt)`. VG
   returns `(VariantCoordinate, HgvsMatchRefAllele, HgvsOriginallyNormalized)`.

5. **`format()` signature**: hgvs_shim uses `format(use_delins_for_inv=False, max_ref_length=None)`.
   VG uses `format(use_compat=False, max_ref_length=settings.HGVS_MAX_REF_ALLELE_LENGTH)`.
   The defaults are equivalent (`False` = don't convert inversions). Only the parameter name differs.

## Architecture Decision

VG's `HGVSConverter` abstract base class keeps its Django-coupled API. The concrete converters
become Django adapters: `PyHGVSConverter` composes `hgvs_shim.PyHGVSConverter` internally;
`BioCommonsHGVSConverter` subclasses `hgvs_shim.BioCommonsHGVSConverter`. Duplicated `*Variant`
classes and exceptions are replaced by direct imports.

---

## Phase 1: Dependency Setup

Add `hgvs_shim` to `requirements.txt` alongside the existing pinned packages — do **not** remove
the git-pinned pyhgvs or biocommons hgvs entries:

```
hgvs_shim==0.2.0
# Keep existing pins:
pyhgvs @ git+https://github.com/SACGF/hgvs@69ff4bdd09088be798708702779a5ba59f1c5b41#egg=pyhgvs
hgvs @ git+https://github.com/biocommons/hgvs@9b32a3077acdb407327c05562fe342d182332c8d#egg=hgvs
```

hgvs_shim declares `hgvs` and (optionally) `pyhgvs` as loose, unpinned dependencies. Pip
will see that your pinned git versions already satisfy those requirements and use them. Your
custom pyhgvs fork and the git-pinned biocommons hgvs will be used everywhere — both by VG
and by hgvs_shim's converters at runtime.

---

## Phase 2: Replace HGVSVariant and Exceptions in hgvs.py

Remove the duplicate class definitions from `genes/hgvs/hgvs.py` and import from `hgvs_shim`:

```python
# Remove these class definitions (~130 lines):
#   HGVSException, HGVSNomenclatureException, HGVSImplementationException, HGVSVariant

# Add at top of genes/hgvs/hgvs.py:
from hgvs_shim import (
    HGVSException,
    HGVSNomenclatureException,
    HGVSImplementationException,
    HGVSVariant,
)
```

### Handle `format()` parameter rename (use_compat → use_delins_for_inv)

The defaults are aligned (`False` in both), so bare `.format()` calls need no changes. Only the
one explicit `use_compat=True` call needs updating:

```python
# classification/models/classification_variant_info_models.py:211
# Before:
c_hgvs_compatible=result.hgvs_variant.format(use_compat=True, max_ref_length=settings.CLASSIFICATION_MAX_REFERENCE_LENGTH),

# After:
c_hgvs_compatible=result.hgvs_variant.format(use_delins_for_inv=True, max_ref_length=settings.CLASSIFICATION_MAX_REFERENCE_LENGTH),
```

---

## Phase 3: Replace PyHGVSVariant; Refactor PyHGVSConverter

### 3.1 `genes/hgvs/pyhgvs/hgvs_converter_pyhgvs.py`

Remove the `PyHGVSVariant` class (~110 lines) and import from hgvs_shim. Refactor
`PyHGVSConverter` to compose `hgvs_shim.PyHGVSConverter` internally.

```python
from hgvs_shim import PyHGVSConverter as _PyHGVSConverter, PyHGVSVariant, TranscriptInfo
from hgvs_shim import HGVSException

class PyHGVSConverter(HGVSConverter):
    """Django-coupled adapter for hgvs_shim.PyHGVSConverter"""

    def __init__(self, genome_build: GenomeBuild, local_resolution=True, clingen_resolution=True):
        super().__init__(genome_build, local_resolution=local_resolution,
                         clingen_resolution=clingen_resolution)
        fasta = genome_build.genome_fasta.fasta
        self._shim = _PyHGVSConverter(fasta_file=fasta, get_transcript=self._load_pyhgvs_transcript)

    def _load_pyhgvs_transcript(self, accession: str):
        """Look up pyhgvs transcript data from the Django DB"""
        from genes.models import TranscriptVersion
        tv = TranscriptVersion.get_transcript_version(self.genome_build, accession)
        return make_transcript(tv.pyhgvs_data)

    @staticmethod
    def _hgvs_name(hgvs_string):
        HGVSConverter._hgvs_string_validation(hgvs_string)
        try:
            return HGVSName(hgvs_string)
        except pyhgvs.InvalidHGVSName as e:
            raise HGVSException(str(e)) from e

    def create_hgvs_variant(self, hgvs_string: str) -> HGVSVariant:
        return self._shim.create_hgvs_variant(hgvs_string)

    def normalize(self, hgvs_variant: PyHGVSVariant) -> HGVSVariant:
        return self._shim.normalize(hgvs_variant)

    def _variant_coordinate_to_g_hgvs(self, vc: VariantCoordinate) -> HGVSVariant:
        chrom, position, ref, alt, _svlen = vc.as_external_explicit(self.genome_build)
        hv = self._shim.variant_coordinate_to_g_hgvs(chrom, position, ref, alt)
        # Fix UCSC chrom → RefSeq contig accession for g.HGVS
        contig = self.genome_build.chrom_contig_mappings[chrom]
        hv._hgvs_name.chrom = contig.refseq_accession
        return hv

    def variant_coordinate_to_c_hgvs(self, vc: VariantCoordinate, transcript_version) -> HGVSVariant:
        chrom, position, ref, alt, _svlen = vc.as_external_explicit(self.genome_build)
        transcript_info = self._tv_to_transcript_info(transcript_version)
        return self._shim.variant_coordinate_to_c_hgvs(chrom, position, ref, alt, transcript_info)

    @staticmethod
    def _tv_to_transcript_info(transcript_version) -> TranscriptInfo:
        """Convert Django TranscriptVersion to hgvs_shim TranscriptInfo"""
        gene_symbol = transcript_version.gene_symbol
        return TranscriptInfo(
            accession=transcript_version.accession,
            strand=transcript_version.strand,
            is_coding=transcript_version.is_coding,
            gene_symbol=str(gene_symbol) if gene_symbol else None,
        )

    def hgvs_to_variant_coordinate_reference_match_and_normalized(
            self, hgvs_string: str, transcript_version=None
    ) -> tuple[VariantCoordinate, HgvsMatchRefAllele, bool]:
        """VG-specific return type wrapping hgvs_shim's simpler hgvs_to_variant_coordinate"""
        pyhgvs_transcript = None
        hgvs_name = self._hgvs_name(hgvs_string)

        if transcript_version:
            pyhgvs_transcript = make_transcript(transcript_version.pyhgvs_data)
            self._validate_in_transcript_range(pyhgvs_transcript, hgvs_name)

        chrom, position, ref, alt = self._shim.hgvs_to_variant_coordinate(hgvs_string)
        contig = self.genome_build.chrom_contig_mappings[chrom]
        chrom = contig.name
        vc = VariantCoordinate.from_explicit_no_svlen(chrom, position, ref, alt)
        matches_reference = self.get_hgvs_match_ref_allele(hgvs_name, pyhgvs_transcript)
        return vc, matches_reference, True

    def c_hgvs_remove_gene_symbol(self, hgvs_string: str) -> str:
        return self._shim.c_hgvs_remove_gene_symbol(hgvs_string)

    def get_transcript_accession(self, hgvs_string: str) -> str:
        return self._shim.get_transcript_accession(hgvs_string)

    def get_hgvs_converter_type(self) -> HGVSConverterType:
        return HGVSConverterType.PYHGVS

    def get_version(self) -> str:
        return self._shim.get_version()

    # --- VG-specific methods ---

    def get_hgvs_match_ref_allele(self, hgvs_name, pyhgvs_transcript=None) -> HgvsMatchRefAllele:
        is_forward_strand = pyhgvs_transcript.tx_position.is_forward_strand if pyhgvs_transcript else True
        ref, _ = hgvs_name.get_ref_alt(is_forward_strand, raw_dup_alleles=True)
        chrom, start, end = hgvs_name.get_raw_coords(pyhgvs_transcript)
        fasta = self.genome_build.genome_fasta.fasta
        genome_ref = get_genomic_sequence(fasta, chrom, start, end)
        return HgvsMatchRefAllele(provided_ref=ref, calculated_ref=genome_ref)

    @staticmethod
    def _validate_in_transcript_range(pyhgvs_transcript, hgvs_name: HGVSName):
        NAMES = {"cdna_start": hgvs_name.cdna_start, "cdna_end": hgvs_name.cdna_end}
        for description, cdna_coord in NAMES.items():
            genomic_coord = pyhgvs_transcript.cdna_to_genomic_coord(cdna_coord)
            within_transcript = (pyhgvs_transcript.tx_position.chrom_start
                                 <= genomic_coord
                                 <= pyhgvs_transcript.tx_position.chrom_stop)
            if not within_transcript:
                raise HGVSException(
                    f"{hgvs_name}: {description} {cdna_coord} resolves outside of transcript")
```

**Note on `hgvs_to_variant_coordinate`**: The shim's call goes through `pyhgvs.parse_hgvs_name()`
using the `get_transcript` callable set in the constructor (a DB lookup). The existing VG code
passed the already-loaded transcript object directly. If this causes a performance issue for
batch use, `hgvs_to_variant_coordinate_reference_match_and_normalized` can call pyhgvs directly
rather than via the shim.

---

## Phase 4: Replace BioCommonsHGVSVariant; Refactor BioCommonsHGVSConverter

### 4.1 `genes/hgvs/biocommons_hgvs/hgvs_converter_biocommons.py`

Remove `BioCommonsHGVSVariant` class and import from hgvs_shim. Subclass
`hgvs_shim.BioCommonsHGVSConverter` to share the biocommons tool setup while keeping VG's
richer `_hgvs_to_g_hgvs` logic (ref fixing, validation, normalization tracking).

Also remove `ParserSingleton` — the latest biocommons HGVS compiles its parser into a Python
class with negligible startup overhead, so the singleton is no longer needed.

```python
from hgvs_shim import BioCommonsHGVSConverter as _BioCommonsHGVSConverterBase
from hgvs_shim import BioCommonsHGVSVariant

class BioCommonsHGVSConverter(_BioCommonsHGVSConverterBase):
    """
    Django-coupled subclass of hgvs_shim.BioCommonsHGVSConverter.
    Adds VariantCoordinate handling, reference matching, normalization tracking,
    and the DjangoTranscriptDataProvider.
    """

    def __init__(self, genome_build: GenomeBuild, local_resolution=True, clingen_resolution=True):
        self.genome_build = genome_build
        self.local_resolution = local_resolution
        self.clingen_resolution = clingen_resolution

        hdp = DjangoTranscriptDataProvider(genome_build)
        if genome_build.name == 'GRCh37':
            assembly_name = genome_build.get_build_with_patch()
        else:
            assembly_name = genome_build.name

        # Calls hgvs_shim base: sets up self.hdp, self.babelfish, self.am, self.ev,
        # self.norm_5p, self.no_validate_mapper, self.no_validate_normalizer
        super().__init__(assembly_name=assembly_name, hdp=hdp)

    def _variant_coordinate_to_sequence_variant(self, vc: VariantCoordinate) -> SequenceVariant:
        chrom, position, ref, alt, _svlen = vc.as_external_explicit(self.genome_build)
        return self.babelfish.vcf_to_g_hgvs(chrom, position, ref, alt)

    def _variant_coordinate_to_g_hgvs(self, vc: VariantCoordinate) -> HGVSVariant:
        var_g = self._variant_coordinate_to_sequence_variant(vc)
        return BioCommonsHGVSVariant(var_g)

    def variant_coordinate_to_c_hgvs(self, vc: VariantCoordinate, transcript_version) -> HGVSVariant:
        """Override to accept TranscriptVersion instead of TranscriptInfo"""
        try:
            var_g = self._variant_coordinate_to_sequence_variant(vc)
            if transcript_version.strand == '-':
                var_g = self.norm_5p.normalize(var_g)
            if transcript_version.is_coding:
                var_c = self.am.g_to_c(var_g, transcript_version.accession)
            else:
                var_c = self.am.g_to_n(var_g, transcript_version.accession)
        except HGVSError as e:
            klass = self._get_exception_class(e)
            raise klass(e) from e

        if gene_symbol := transcript_version.gene_symbol:
            var_c.gene = gene_symbol.symbol
        return BioCommonsHGVSVariant(var_c)

    def hgvs_to_variant_coordinate_reference_match_and_normalized(
            self, hgvs_string: str, transcript_version=None
    ) -> tuple[VariantCoordinate, HgvsMatchRefAllele, HgvsOriginallyNormalized]:
        # Keep existing VG implementation (calls self._hgvs_to_g_hgvs which stays in VG)
        try:
            var_g, matches_reference, originally_normalized = self._hgvs_to_g_hgvs(hgvs_string)
            try:
                (chrom, position, ref, alt, _typ) = self.babelfish.hgvs_to_vcf(var_g)
                if alt == '.':
                    alt = ref
            except HGVSDataNotAvailableError:
                raise Contig.ContigNotInBuildError()
        except HGVSError as hgvs_error:
            klass = self._get_exception_class(hgvs_error)
            raise klass(hgvs_error) from hgvs_error

        vc = VariantCoordinate.from_explicit_no_svlen(chrom, position, ref=ref, alt=alt)
        return vc.as_internal_symbolic(self.genome_build), matches_reference, originally_normalized

    def get_hgvs_converter_type(self) -> HGVSConverterType:
        return HGVSConverterType.BIOCOMMONS_HGVS

    def description(self, describe_fallback=True) -> str:
        hgvs_converter_type = self.get_hgvs_converter_type()
        version = self.get_version()
        desc = f"{hgvs_converter_type.name} {version}"
        if describe_fallback and self.clingen_resolution:
            desc += " (ClinGen fallback)"
        return desc

    # Keep VG-specific private methods unchanged:
    # _fix_ref(), _hgvs_to_g_hgvs() (complex version with ref fixing + validation),
    # HgvsMatchTranscriptAndGenomeRefAllele
```

### 4.2 Keep `HgvsMatchTranscriptAndGenomeRefAllele` in VG

VG-specific class with strand-aware ref comparison. Stays in `hgvs_converter_biocommons.py`.

---

## Phase 5: ComboCheckerHGVSConverter — No Change

VG's version uses `VariantCoordinate` in method signatures; hgvs_shim's uses `(chrom, pos, ref, alt)`.
VG's version is correct for VG's internal API. No changes needed.

---

## Phase 6: Update HGVSConverterType

hgvs_shim 0.2.0 uses `IntEnum`, so VG can define its own extended enum safely. Keep VG's full
definition in `hgvs_converter.py` — it's small and avoids any import-order issues:

```python
# genes/hgvs/hgvs_converter.py
# Keep existing VG HGVSConverterType as-is — values 1,2,3 match hgvs_shim's IntEnum values.
# CLINGEN_ALLELE_REGISTRY = 4 is VG-specific.
class HGVSConverterType(Enum):
    PYHGVS = 1
    BIOCOMMONS_HGVS = 2
    COMBO = 3
    CLINGEN_ALLELE_REGISTRY = 4
    ...
```

---

## Phase 7: Update `genes/hgvs/__init__.py`

```python
from hgvs_shim import HGVSException, HGVSNomenclatureException, HGVSImplementationException
from hgvs_shim import HGVSVariant
from .hgvs import *       # CHGVS, PHGVS, CHGVSDiff, chgvs_diff_description
from .hgvs_matcher import *  # HGVSMatcher, VariantResolvingError, etc.
```

---

## Phase 8: data_provider.py — No Change

`genes/hgvs/biocommons_hgvs/data_provider.py` is fully Django-specific with no hgvs_shim
equivalent. Keep as-is.

---

## Phase 9: Testing

Run after each phase:

```bash
python3 manage.py test --keepdb genes.tests.test_hgvs
python3 manage.py test --keepdb genes.tests.test_clean_hgvs
python3 manage.py test --keepdb genes
python3 manage.py test --keepdb classification
```

Key things to verify:
- g.HGVS contig accession mapping still works after PyHGVSConverter refactor
- BioCommonsHGVSConverter `_hgvs_to_g_hgvs` validation logic still works after subclassing
- Inversion formatting unchanged (both old and new default to `False`)

---

## File-by-file Change Summary

| File | Change |
|---|---|
| `requirements.txt` | Add `hgvs_shim==0.2.0` (keep existing pinned pyhgvs and hgvs) |
| `genes/hgvs/hgvs.py` | Remove `HGVSVariant`, `HGVSException*`; import from hgvs_shim |
| `genes/hgvs/__init__.py` | Update re-exports to import from hgvs_shim |
| `genes/hgvs/hgvs_converter.py` | No change to `HGVSConverterType`; minor import updates |
| `genes/hgvs/pyhgvs/hgvs_converter_pyhgvs.py` | Remove `PyHGVSVariant`; refactor `PyHGVSConverter` to compose shim |
| `genes/hgvs/biocommons_hgvs/hgvs_converter_biocommons.py` | Remove `BioCommonsHGVSVariant`, `ParserSingleton`; subclass shim |
| `genes/hgvs/hgvs_converter_combo.py` | No change |
| `genes/hgvs/biocommons_hgvs/data_provider.py` | No change |
| `genes/hgvs/hgvs_matcher.py` | No change |
| `classification/models/classification_variant_info_models.py` | `use_compat=True` → `use_delins_for_inv=True` (one call site) |

---

## Migration Checklist

- [ ] Add `hgvs_shim==0.2.0` to `requirements.txt` (keep existing pinned pyhgvs and hgvs)
- [ ] Remove `HGVSException*` and `HGVSVariant` from `hgvs.py`, import from hgvs_shim
- [ ] Update `genes/hgvs/__init__.py` re-exports
- [ ] Remove `PyHGVSVariant` from `hgvs_converter_pyhgvs.py`, import from hgvs_shim
- [ ] Refactor `PyHGVSConverter` to compose `hgvs_shim.PyHGVSConverter`
- [ ] Remove `BioCommonsHGVSVariant` and `ParserSingleton` from `hgvs_converter_biocommons.py`
- [ ] Import `BioCommonsHGVSVariant` from hgvs_shim
- [ ] Refactor `BioCommonsHGVSConverter` to subclass `hgvs_shim.BioCommonsHGVSConverter`
- [ ] Fix `use_compat=True` → `use_delins_for_inv=True` in `classification_variant_info_models.py`
- [ ] Run test suite: `python3 manage.py test --keepdb genes classification`

---

## Implementation Order

1. `requirements.txt` (no code changes, verify install works)
2. Exceptions — safest, no interface change
3. `HGVSVariant` — abstract, no runtime impact
4. `PyHGVSVariant` / `BioCommonsHGVSVariant` — simple wrappers
5. `PyHGVSConverter` — test carefully (chrom mapping, ref matching)
6. `BioCommonsHGVSConverter` — largest change, test thoroughly

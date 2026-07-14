# Plan: Speed up SV annotation by moving conservation bigWig scoring out of VEP into pyBigWig

Tracking issue: [#1657](https://github.com/SACGF/variantgrid/issues/1657)

## Background & finding

Structural-variant (SV) annotation runs can take hours / never finish. The suspicion was
the gnomAD-SV overlap. **Benchmarking disproved that** — the bottleneck is the four
conservation **bigWig** `--custom` annotations (phastCons/phyloP), whose cost scales with
**SV span length**.

### What the SV VEP command actually runs

`get_vep_command(..., pipeline_type=STRUCTURAL_VARIANT)` (`annotation/vep_annotation.py`)
emits, per SV, these `--custom` overlaps (GRCh37 RefSeq shown):

- `gnomAD_SV`, `gnomAD_SV_name` — VCF overlap, `overlap_cutoff=80,same_type=1`
- `phastCons100way`, `phastCons46way`, `phyloP100way`, `phyloP46way` —
  `type=overlap,summary_stats=max,format=bigwig` (5–9 GB `.bw` files each)
- `REPEAT_MASKER` — bed overlap

The 4 conservation tracks are configured (via `ColumnVEPField` / `vep_columns.filter_for`)
to apply to **both** STANDARD and STRUCTURAL_VARIANT pipelines. On 1 bp small variants they
are cheap; on multi-Mb SVs `summary_stats=max` forces VEP's Perl bigWig reader to walk
conservation values across the entire SV span, ×4 tracks.

### Benchmark evidence (GRCh37 / RefSeq / VEP 116, real prod data files)

Query A — 1,486 SVs sampled from gnomAD (mostly small, some large):

| VEP configuration | time |
|---|---|
| baseline (transcript consequence, no `--custom`) | 10.4 s |
| + gnomAD_SV overlap (2 customs) | 13.3 s (**+2.9 s**) |
| + RepeatMasker (bed) | 11.6 s (**+1.2 s**) |
| + 4 conservation bigwigs | **timeout >290 s** |
| full SV command (all customs) | **timeout >120 s** |

Span dependence (single bigwig track):

| query | span | time |
|---|---|---|
| 559 tiny INS | ~50 bp | 5.6 s |
| 200 large SVs | 0.5–4.5 Mb | **timeout >200 s** |

Conclusion: gnomAD overlap and RepeatMasker are trivial; the conservation bigwigs are the
bottleneck and the cost is **O(SV span)**. vcfanno does **not** apply here — it can't read
bigWig, and the slow thing isn't the gnomAD VCF overlap.

## End result / target

Drop the 4 conservation `--custom` bigwig args from the **SV** VEP command, and compute the
same 4 conservation `_max` columns in a separate **pyBigWig** stage. Expected: SV VEP goes
from *never-finishes* to ~baseline (~8–11 s even on large SVs); conservation added in well
under a second for a normal batch; values within tolerance of VEP.

Measured PoC vs VEP (200 large SVs, 4 tracks): VEP **>290 s (never finished)** →
pyBigWig `aligned` **10.2 s, bit-exact** (err 0.000 vs per-base ground truth).

## Key implementation facts discovered

### 1. VEP's summary window convention

VEP takes `summary_stats=max` over the SV's reference footprint =
**1-based inclusive `[POS, max(END, POS+|SVLEN|)]`** = 0-based half-open
`[POS-1, max(END, POS+|SVLEN|))`:

- DEL/DUP: `END` (== `POS+|SVLEN|`)
- INS: `END` is only `POS+1`, so the footprint is `POS+|SVLEN|`

Validated: with this window, pyBigWig `exact` matched VEP within 0.05 on 100/100/99.7/100 %
of moderate SVs (phastCons100/phastCons46/phyloP100/phyloP46). The lone miss is a 0.40
discrepancy on the VEP side.

### 2. Fast + exact summary: aligned-window decomposition

Naive `pyBigWig.stats(exact=False)` (zoom summaries) is ~500× faster than per-base but
**overshoots MAX** at span edges (phyloP off by up to 1.27), because zoom bins straddle the
query boundary and pull in out-of-span values.

The `aligned` decomposition is fast **and** bit-exact:

- exact per-base on the two ragged ends (width `E`)
- interior read from zoom summaries on a fine grid (bin width `G`), one `stats(nBins=...)` call
- **safety invariant:** `E > G` (ends wider than the interior's zoom reduction), so any
  straddle from an interior bin lands inside an exact-scanned end and is capped correctly

Chosen `G = 16384`, `E = 2*G`. Spans `<= 2*E` just do a full exact scan.

Speed / accuracy across modes (200 large SVs, 0.5–4.5 Mb, 4 tracks):

| mode | large-SV time | vs exact ground truth |
|---|---|---|
| VEP `--custom` bigwig | >290 s (never finished) | — |
| pyBigWig `exact` (per-base) | 41.3 s | ground truth |
| pyBigWig **`aligned`** | **10.2 s** | **100 % exact (err 0.000)** |
| pyBigWig `hybrid` (128k exact ends) | 3.8 s | phyloP 98 % (max 0.56) |
| pyBigWig `zoom` | 0.09 s | overshoots (phyloP up to 1.27) |

On normal SVs (≤30 kb) every pyBigWig mode is ~1 ms/SV vs VEP's 293 ms/SV.

`aligned` also produces `average` from the same decomposition (coverage-weighted combine of
interior zoom-mean bins + exact ends); validated vs exact mean, max error 0.012 on large SVs.
Production only needs `max`, but the method covers both.

## Code sketch — the pyBigWig summary (from PoC `claude/sv_pybigwig_poc.py`, mode `aligned`)

```python
import pyBigWig

def span_max_aligned(bw, chrom, s, e):
    """MAX over 0-based half-open [s, e) — bit-exact, fast on huge spans.
    exact per-base on ragged ends (width E), interior from zoom summaries on a fine
    aligned grid (bin width G). Invariant E > G caps any zoom straddle in the exact ends."""
    G = 16384          # interior zoom bin width -> zoom reduction used is <= G
    E = 2 * G          # exact-end width must exceed that reduction
    if e - s <= 2 * E:
        return bw.stats(chrom, s, e, type="max", exact=True)[0]
    lo, hi = s + E, e - E
    vals = [
        bw.stats(chrom, s, lo, type="max", exact=True)[0],
        bw.stats(chrom, hi, e, type="max", exact=True)[0],
    ]
    nb = max(1, (hi - lo) // G)
    vals.extend(v for v in bw.stats(chrom, lo, hi, type="max", exact=False, nBins=nb)
                if v is not None)
    vals = [v for v in vals if v is not None]
    return max(vals) if vals else None


# per-variant window (VEP-equivalent footprint)
def sv_window(pos, end, svlen):        # svlen may be None
    s = pos - 1                        # 0-based inclusive start
    e = max(end, pos + abs(svlen or 0))  # 0-based exclusive end (INS uses POS+|SVLEN|)
    return s, e

# chrom naming differs per file: ck = ("chr" if "chr1" in bw.chroms() else "") + chrom
# clamp: s>=0, e<=chromlen, e>s
```

The PoC also implements `exact` / `zoom` / `hybrid` modes and the mean/average variant for
comparison — see `claude/sv_pybigwig_poc.py`.

## Integration sketch (to design later)

1. **Stop emitting the 4 conservation bigwig customs on the SV pipeline.** They apply to SV
   only because their `ColumnVEPField` rows are configured for STRUCTURAL_VARIANT; restrict
   them to STANDARD so `get_vep_command(..., STRUCTURAL_VARIANT)` drops them
   (`vep_columns.filter_for(pipeline_type=...)` already gates this). Small-variant STANDARD
   behaviour unchanged.

2. **Add a post-VEP pyBigWig stage** in `dump_and_annotate_variants`
   (`annotation/tasks/annotate_variants.py`), alongside the existing AnnotSV stage, for
   STRUCTURAL_VARIANT runs: read the run's dumped SV VCF, compute the 4 conservation `_max`
   columns with `span_max_aligned`, and merge them into the annotated VCF / import path so
   the columns land in the same `_max` fields VEP would have written. Parallelizable with VEP.

3. **Settings / config** — bigWig paths already live in `VEPConfig`
   (`phastcons100way`, `phastcons46way`, `phylop100way`, `phylop46way`). Gate the new stage
   behind a setting; add `pyBigWig` as a dependency.

4. **Test** — diff pyBigWig output vs a VEP-annotated small SV VCF (fixtures already exist,
   e.g. `annotation/tests/test_data/test_columns_version*_grch3*_sv.vep_annotated.vcf`),
   asserting the 4 `_max` columns agree within tolerance (~0.05).

## Caveats

- Validated on GRCh37 / RefSeq / VEP 116. Re-check the window convention and the
  `E`/`G` straddle-safety on GRCh38 (gnomAD v4 SV file is larger, but that's the cheap
  overlap, not the bigwig issue).
- `aligned` is bit-exact; `hybrid` trades a ~0.5 phyloP tolerance for ~3× more speed.
  Recommend shipping `aligned`.
- Conservation `max` over a multi-Mb SV is saturated/near-meaningless anyway (phastCons ≈ 1.0),
  so an alternative worth considering later is dropping these columns for SVs entirely rather
  than recomputing them — but this plan preserves current output.

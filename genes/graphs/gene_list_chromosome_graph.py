import logging
from collections import defaultdict

from matplotlib.patches import Rectangle

from genes.models import GeneList
from library.graphs.chromosomes_graph import read_cytoband, centromere_mid
from library.utils import sha256sum_str
from snpdb.graphs.graphcache import CacheableGraph
from snpdb.models.models_enums import AssemblyMoleculeType
from snpdb.models.models_genome import GenomeBuild

CHROMOSOME_WIDTH = 0.07      # data units – fraction of spacing slot
SCALE = 1e6                  # 1 data unit = 1 Mb (y axis)
GENE_LABEL_FONT_SIZE = 7
GENE_LABEL_FONT_WEIGHT = 'bold'
GENE_LABEL_GAP = 0.005       # gap between chromosome edge and gene label text (data units)
CHROM_LABEL_FONT_SIZE = 9
CHROM_LABEL_PADDING = 3.0   # Mb below the lowest chromosome bottom in each row
NUM_ROWS = 3

# Spacing calculation constants
_FIG_WIDTH_IN = 14.0
_FIG_USABLE_FRAC = 0.98     # left=0.01, right=0.99
_CHAR_WIDTH_FACTOR = 0.55   # character width as fraction of font size (in points)
_MIN_INTER_LABEL_GAP_IN = 0.10  # minimum gap between adjacent gene labels (inches)
_MIN_CHROM_SPACING = 0.18   # ~39% fill ratio (CW/spacing); formula expands for long names

# Simplified banding colours
BAND_COLORS = {
    'gneg':    (.95, .95, .95),
    'gpos25':  (.75, .75, .75),
    'gpos50':  (.50, .50, .50),
    'gpos75':  (.30, .30, .30),
    'gpos100': (.10, .10, .10),
    'acen':    (.78, .25, .25),
    'gvar':    (.70, .70, .70),
    'stalk':   (.85, .85, .85),
}

MINUS_COLOR = 'blue'   # '-' strand → left
PLUS_COLOR  = 'red'    # '+' strand → right


class GeneListChromosomeGraph(CacheableGraph):
    """
    Karyotype-style chromosome graph for a gene list.

    Layout
    ------
    * Chromosomes are arranged in three rows (1-8, 9-16, 17-22+X+Y).
    * Each chromosome is drawn as a vertical bar with simplified G-banding.
    * Centromeres are aligned to a dashed horizontal baseline on each row.
    * p arm extends upward (shorter arm, at the top).
    * q arm extends downward.
    * No numeric axis – each chromosome is labelled below its bar.
    * '-' strand genes: blue tick + label to the *left*.
    * '+' strand genes: red  tick + label to the *right*.
    """

    def __init__(self, gene_list_id):
        super().__init__()
        self.gene_list_id = gene_list_id
        self.figsize = (14, 12)

    def get_params_hash(self):
        return sha256sum_str(str(self.gene_list_id))

    # ------------------------------------------------------------------ #
    # CacheableGraph interface                                             #
    # ------------------------------------------------------------------ #

    def plot(self, ax):
        """Not used – plot_figure() drives everything."""

    @staticmethod
    def _compute_spacing(rows, gene_data):
        """
        Compute a single uniform inter-chromosome spacing for all rows so that
        adjacent gene labels just don't touch anywhere, solved analytically.
        Uses the longest row as n (most constraining for the formula).
        """
        n = max(len(row) for row in rows)
        W = _FIG_WIDTH_IN * _FIG_USABLE_FRAC
        char_w_in = GENE_LABEL_FONT_SIZE * _CHAR_WIDTH_FACTOR / 72

        # Worst-case adjacent pair across ALL rows
        max_pair_chars = 0
        for row_chroms in rows:
            for i in range(len(row_chroms) - 1):
                right_max = max(
                    (len(nm) for _, nm, s in gene_data.get(row_chroms[i], []) if s != '-'),
                    default=0,
                )
                left_max = max(
                    (len(nm) for _, nm, s in gene_data.get(row_chroms[i + 1], []) if s == '-'),
                    default=0,
                )
                max_pair_chars = max(max_pair_chars, right_max + left_max)

        label_gap_in = max_pair_chars * char_w_in + _MIN_INTER_LABEL_GAP_IN
        factor = n * label_gap_in / W
        if factor >= 0.95:
            return CHROMOSOME_WIDTH + label_gap_in  # degenerate fallback
        spacing = CHROMOSOME_WIDTH / (1 - factor)
        return max(spacing, _MIN_CHROM_SPACING)

    @staticmethod
    def _compute_xlim(rows, gene_data, spacing):
        """
        Compute a shared xlim for all rows so chromosomes are horizontally aligned.
        Based on the widest row (most chromosomes) and worst outer labels.
        """
        n_max = max(len(row) for row in rows)
        # Approximate: 1 char ≈ font_size * char_factor / 72 inches; convert to data units
        char_w_data = (GENE_LABEL_FONT_SIZE * _CHAR_WIDTH_FACTOR / 72) * \
                      (n_max * spacing) / (_FIG_WIDTH_IN * _FIG_USABLE_FRAC)
        _edge = CHROMOSOME_WIDTH / 2 + GENE_LABEL_GAP

        left_max = _edge
        right_max = _edge
        for row_chroms in rows:
            left = max(
                (len(nm) * char_w_data + _edge
                 for _, nm, s in gene_data.get(row_chroms[0], []) if s == '-'),
                default=_edge,
            )
            right = max(
                (len(nm) * char_w_data + _edge
                 for _, nm, s in gene_data.get(row_chroms[-1], []) if s != '-'),
                default=_edge,
            )
            left_max = max(left_max, left)
            right_max = max(right_max, right)

        return -left_max, (n_max - 1) * spacing + right_max

    @staticmethod
    def _build_rows(genome_build):
        """
        Divide the genome build's chromosomes into NUM_ROWS rows.
        Preserves the build's natural ordering; remainder goes on the last row.
        """
        chroms = list(
            genome_build.standard_contigs
            .filter(molecule_type=AssemblyMoleculeType.CHROMOSOME)
            .values_list("name", flat=True)
        )
        n = len(chroms)
        base = n // NUM_ROWS
        rows = [chroms[i * base:(i + 1) * base] for i in range(NUM_ROWS - 1)]
        rows.append(chroms[(NUM_ROWS - 1) * base:])  # last row gets remainder
        return [row for row in rows if row]

    def plot_figure(self, figure):
        genome_build = GenomeBuild.grch38()
        cytoband_data = read_cytoband(genome_build.settings["cytoband"])
        gene_data = self._collect_gene_data(genome_build)
        rows = self._build_rows(genome_build)

        # Work out proportional row heights (based on tallest chromosome in row)
        row_heights = []
        for row_chroms in rows:
            max_len = 0
            for ch in row_chroms:
                bands = cytoband_data.get(f"chr{ch}", [])
                if bands:
                    max_len = max(max_len, bands[-1][1] / SCALE)  # end in Mb
            row_heights.append(max(max_len, 10))  # at least 10 Mb

        spacing = self._compute_spacing(rows, gene_data)
        xlim = self._compute_xlim(rows, gene_data, spacing)

        gs = figure.add_gridspec(len(rows), 1,
                                  height_ratios=row_heights,
                                  hspace=0.15)

        for row_idx, row_chroms in enumerate(rows):
            ax = figure.add_subplot(gs[row_idx])
            self._draw_row(ax, row_chroms, cytoband_data, gene_data, spacing, xlim)

    # ------------------------------------------------------------------ #
    # Data loading                                                         #
    # ------------------------------------------------------------------ #

    def _collect_gene_data(self, genome_build):
        """Return dict: bare_chrom ('1','X'…) → list of (pos_bp, name, strand)."""
        gene_list = GeneList.objects.get(pk=self.gene_list_id)
        chrom_genes = defaultdict(list)

        for glgs in gene_list.genelistgenesymbol_set.filter(
                gene_symbol__isnull=False).select_related('gene_symbol'):
            gs = glgs.gene_symbol
            gv = gs.latest_gene_version(genome_build)
            if gv is None:
                continue
            try:
                chrom = gv.chrom        # '1', '17', 'X' etc. (no 'chr' prefix)
                midpoint = (gv.start + gv.end) / 2
                strand = gv.strand
            except Exception:
                logging.warning("Could not get coordinates for %s", gs)
                continue
            chrom_genes[chrom].append((midpoint, str(gs), strand))

        return chrom_genes

    # ------------------------------------------------------------------ #
    # Drawing                                                              #
    # ------------------------------------------------------------------ #

    def _draw_row(self, ax, row_chroms, cytoband_data, gene_data, spacing, xlim):
        max_p_mb = 0   # max p-arm length across this row
        max_q_mb = 0   # max q-arm length across this row

        # First pass: gather extents so we can set ylim correctly
        row_info = {}   # chrom_num → (cen_mid, bands)
        for chrom_num in row_chroms:
            bands = cytoband_data.get(f"chr{chrom_num}", [])
            if not bands:
                continue
            cen_mid = centromere_mid(bands)
            p_mb = cen_mid / SCALE
            q_mb = (bands[-1][1] - cen_mid) / SCALE
            max_p_mb = max(max_p_mb, p_mb)
            max_q_mb = max(max_q_mb, q_mb)
            row_info[chrom_num] = (cen_mid, bands)

        # All chromosome labels share the same y: bottom of longest q-arm + padding
        row_label_y = -(max_q_mb + CHROM_LABEL_PADDING)

        # Second pass: draw
        for x_idx, chrom_num in enumerate(row_chroms):
            if chrom_num not in row_info:
                continue
            cen_mid, bands = row_info[chrom_num]
            x = x_idx * spacing

            # --- Bands ---
            for b_start, b_end, _name, stain in bands:
                y0 = -(b_start - cen_mid) / SCALE   # inverted: p arm upward
                y1 = -(b_end   - cen_mid) / SCALE
                yb, yt = min(y0, y1), max(y0, y1)
                color = BAND_COLORS.get(stain, (.6, .6, .6))
                rect = Rectangle(
                    (x - CHROMOSOME_WIDTH / 2, yb),
                    CHROMOSOME_WIDTH, yt - yb,
                    facecolor=color, edgecolor='none', linewidth=0,
                )
                ax.add_patch(rect)

            # Chromosome outline
            chrom_top    = -(bands[0][0]  - cen_mid) / SCALE
            chrom_bottom = -(bands[-1][1] - cen_mid) / SCALE
            outline = Rectangle(
                (x - CHROMOSOME_WIDTH / 2, chrom_bottom),
                CHROMOSOME_WIDTH, chrom_top - chrom_bottom,
                facecolor='none', edgecolor='black', linewidth=0.4,
            )
            ax.add_patch(outline)

            # --- Chromosome label (all aligned to same y in this row) ---
            ax.text(x, row_label_y, chrom_num,
                    ha='center', va='top',
                    fontsize=CHROM_LABEL_FONT_SIZE, color='black')

            # --- Gene labels (no tick line, placed directly beside chromosome) ---
            for pos_bp, gene_name, strand in gene_data.get(chrom_num, []):
                gene_y = -(pos_bp - cen_mid) / SCALE   # same inversion

                if strand == '-':
                    ax.text(x - CHROMOSOME_WIDTH / 2 - GENE_LABEL_GAP, gene_y,
                            gene_name, ha='right', va='center',
                            fontsize=GENE_LABEL_FONT_SIZE,
                            fontweight=GENE_LABEL_FONT_WEIGHT,
                            color=MINUS_COLOR, clip_on=False)
                else:
                    ax.text(x + CHROMOSOME_WIDTH / 2 + GENE_LABEL_GAP, gene_y,
                            gene_name, ha='left', va='center',
                            fontsize=GENE_LABEL_FONT_SIZE,
                            fontweight=GENE_LABEL_FONT_WEIGHT,
                            color=PLUS_COLOR, clip_on=False)

        # --- Centromere dashed baseline ---
        ax.axhline(0, color='#888888', linestyle='--', linewidth=0.7, alpha=0.8)

        # --- Axis limits and styling ---
        ax.set_xlim(*xlim)
        ax.set_ylim(row_label_y - 2, max_p_mb + 1)
        ax.axis('off')

    def figure(self, figure):
        figure.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)

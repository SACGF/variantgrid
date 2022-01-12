from shlex import shlex

from library import file_utils, genomics

BED_DETAILS = 'bedDetail'
BED_FILE = 'bedFile'
REQUIRED_FIELDS = ['chrom', 'chrom_start', 'chrom_end']
OPTIONAL_FIELDS = ['name', 'score', 'strand']
THICK_ATTRIBUTES = ['thick_start', 'thick_end']
RGB_ATTRIBUTE = ['item_rgb']
BLOCK_ATTRIBUTES = ['block_count', 'block_sizes', 'block_starts']
OPTIONAL_ATTRIBUTES = THICK_ATTRIBUTES + RGB_ATTRIBUTE + BLOCK_ATTRIBUTES
ALL_OPTIONAL = OPTIONAL_FIELDS + OPTIONAL_ATTRIBUTES

TYPE_CONVERSIONS = {'chrom_start': int,
                    'chrom_end': int,
                    'score': float,
                    'thick_start': int,
                    'thick_end': int,
                    'block_count': int,
                    'block_sizes': str,
                    'block_starts': str}


class BedFileReader:
    """
    Reads bed file
    Returns - iterator of HTSeq Genomic Features
    Attribute fields (column 7-12: 'thick_start', 'thick_end', 'block_count', 'block_sizes', 'block_starts')
        are saved as dict in feature attribute 'attr'. e.g. feature.attr['thick_start']
    """

    def __init__(self, filename, **kwargs):
        """ @param want_chr: set to True/False to FORCE add/removal of "chr" (default = leave as is) """
        self.expected_genome = kwargs.get("expected_genome")
        self.source = file_utils.file_or_filename(filename, "rt")
        self.delimiter = kwargs.get("delimiter", '\t')
        self.want_chr = kwargs.get("want_chr")
        self.bed_type = BED_FILE
        self.line_buffer = []
        self.track = {}
        self.read_header()

    def track_to_dict(self, track_line):
        data = {}
        try:
            lexer = shlex(track_line, posix=True)
            lexer.whitespace += ' '
            lexer.wordchars += '='
            lines = list(lexer)
            if lines[0] != 'track':
                raise ValueError(f"Not a track line: '{track_line}'")
            for word in lines[1:]:
                k, v = word.split("=", 1)
                data[k] = v
        except:
            pass
        return data

    def get_genome_build_name(self):
        POTENTIAL_GENOME_BUILD_NAMES = ["db", "genome", "genome_build", "genomeBuild"]
        if self.track:
            for name in POTENTIAL_GENOME_BUILD_NAMES:
                val = self.track.get(name)
                if val:
                    return val
        return False

    def read_header(self):
        """ Read until we find a non-header line. Save this as line_buffer """
        while True:
            line = next(self.source)
            if line.startswith('track'):
                self.track = self.track_to_dict(line)
                bed_type = self.track.get("type")
                if bed_type == BED_DETAILS:
                    self.bed_type = BED_DETAILS
            elif line.startswith('#') or line.startswith('browser'):
                pass
            else:
                self.line_buffer.append(line)
                break

    def __iter__(self):
        return self

    def __next__(self):  # @ReservedAssignment
        # Read until we get a non-comment line
        while True:
            if self.line_buffer:
                line = self.line_buffer.pop(0)
            else:
                line = next(self.source)
            if not line.startswith("#"):
                break

        columns = line.rstrip().split(self.delimiter)
        details_columns = []
        if self.bed_type == BED_DETAILS:
            num_columns = len(columns)
            if num_columns != 14:
                raise ValueError(f"bedDetails line should have 14 columns, was {num_columns}: '{line}'")

            # Last 2 fields are details - strip them
            columns = columns[:-2]
            details_columns = columns[-2:]

        data = dict(list(zip(REQUIRED_FIELDS, columns[:len(REQUIRED_FIELDS)])))
        # Only put in optional fields when provided
        data.update(zip(ALL_OPTIONAL, columns[len(REQUIRED_FIELDS):]))
        # Convert to appropriate types
        for f in REQUIRED_FIELDS + ALL_OPTIONAL:
            if f in TYPE_CONVERSIONS and data.get(f):
                convert = TYPE_CONVERSIONS[f]
                data[f] = convert(data[f])

        # Make sure they aren't negative (some dodgy .BED files around)
        start = max(0, data["chrom_start"])
        end = max(0, data["chrom_end"])
        chrom = data["chrom"]
        if self.want_chr is not None:  # format chrom
            chrom = genomics.format_chrom(chrom, self.want_chr)
        import HTSeq
        iv = HTSeq.GenomicInterval(chrom, start, end)
        name = data.get("name", "")
        score = data.get("score", 0.0)  # Using "." here causes intron to not be displayed correctly
        iv.strand = data.get("strand", '.')

        feature = HTSeq.GenomicFeature(name, "from_bed", iv)
        feature.score = score
        feature.attr = {"score": score}  # So attr and .score are always the same (may be 0.0 if not found)

        # Copy set optional fields to GenomicFeature.attr
        for field in OPTIONAL_ATTRIBUTES:
            value = data.get(field)
            if value:
                feature.attr[field] = value

        if details_columns:
            feature.attr["ID"] = details_columns[0]
            feature.attr["description"] = details_columns[1]

        return feature

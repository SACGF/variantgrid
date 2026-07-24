from rest_framework import serializers

from genes.models import GeneInfo, GeneListCategory, GeneList, Gene, Transcript, GeneListGeneSymbol, \
    GeneAnnotationRelease, SampleGeneList, ActiveSampleGeneList, GeneSymbol, TranscriptVersion, GeneVersion, HGNC, \
    GeneCoverageCollection, GeneCoverageCanonicalTranscript
from snpdb.models import Company, Contig
from snpdb.serializers import GenomeBuildSerializer


class GeneSymbolSerializer(serializers.ModelSerializer):

    class Meta:
        model = GeneSymbol
        fields = ('symbol', )


class GeneSerializer(serializers.ModelSerializer):

    class Meta:
        model = Gene
        fields = ('identifier', )


class GeneVersionSerializer(serializers.ModelSerializer):
    gene = GeneSerializer()
    gene_symbol = GeneSymbolSerializer()
    genome_build = GenomeBuildSerializer()

    class Meta:
        model = GeneVersion
        fields = ('gene', 'version', 'gene_symbol', 'genome_build')


class TranscriptSerializer(serializers.ModelSerializer):

    class Meta:
        model = Transcript
        fields = ('identifier', )


class TranscriptVersionSerializer(serializers.ModelSerializer):
    transcript = TranscriptSerializer()
    genome_build = GenomeBuildSerializer()
    gene_version = GeneVersionSerializer()

    class Meta:
        model = TranscriptVersion
        fields = ('transcript', 'version', 'genome_build', 'gene_version')


class HGNCSerializer(serializers.ModelSerializer):
    gene_symbol = GeneSymbolSerializer()
    status = serializers.SerializerMethodField()

    class Meta:
        model = HGNC
        fields = ('hgnc_id', 'gene_symbol', 'approved_name', 'status')

    def get_status(self, obj: HGNC):
        return obj.get_status_display()


class GeneListCategorySerializer(serializers.ModelSerializer):

    class Meta:
        model = GeneListCategory
        fields = ('name', 'public', 'company', 'icon_css_class', 'description')


class GeneListGeneSymbolSerializer(serializers.ModelSerializer):
    class Meta:
        model = GeneListGeneSymbol
        fields = ('original_name', 'gene_symbol', 'gene_symbol_alias', 'modification_info')


class GeneListSerializer(serializers.ModelSerializer):
    category = GeneListCategorySerializer(allow_null=True)
    user = serializers.StringRelatedField()
    genelistgenesymbol_set = GeneListGeneSymbolSerializer(many=True)
    can_write = serializers.SerializerMethodField()
    absolute_url = serializers.URLField(source='get_absolute_url', read_only=True)

    class Meta:
        model = GeneList
        fields = ('pk', 'category', 'name', 'user', 'import_status', 'genelistgenesymbol_set', 'can_write', 'absolute_url')

    def get_can_write(self, obj: GeneList) -> bool:
        request = self.context.get('request')
        if request is not None:
            user = request.user
        else:
            user = self.context.get('user')
        if user is None:
            return False

        if obj.id:  # Can't check pk as may be FakeGeneList object
            can_write = obj.can_write(user)
        else:
            # Not yet saved, so no permissions etc
            can_write = user == obj.user

        if can_write and obj.category and obj.category.company:
            company = Company.get_our_company()
            if obj.category.company != company:
                can_write = False

        return can_write

    def create(self, validated_data):
        genelistgenesymbol_set_data = validated_data.pop('genelistgenesymbol_set', [])
        validated_data.pop('category', None)  # System-specific FK, skip on import
        user = self.context.get('user')
        gene_list = GeneList.objects.create(user=user, **validated_data)
        for gs_data in genelistgenesymbol_set_data:
            gs_data.pop('gene_symbol_alias', None)  # System-specific integer PK, skip on import
            GeneListGeneSymbol.objects.create(gene_list=gene_list, **gs_data)
        return gene_list

    def get_fields(self):
        # Code from https://timmytango.com/notes/excluding-fields-in-drf-serializers/
        fields = super().get_fields()

        exclude_fields = self.context.get('exclude_fields', [])
        for field in exclude_fields:
            # providing a default prevents a KeyError
            # if the field does not exist
            fields.pop(field, None)

        return fields


class GeneAnnotationReleaseSerializer(serializers.ModelSerializer):
    __str__ = serializers.SerializerMethodField()

    class Meta:
        model = GeneAnnotationRelease
        fields = ("version", "annotation_consortium", "genome_build", "__str__")

    def get___str__(self, obj) -> str:
        return str(obj)


class GeneContigSerializer(serializers.ModelSerializer):
    """ Includes the pk, so callers can filter variant queries on locus__contig_id """

    class Meta:
        model = Contig
        fields = ('id', 'name', 'refseq_accession')


class GenomicCoordinatesSerializer(serializers.Serializer):
    """ Serializes a GenomicCoordinates dataclass """
    chrom = serializers.CharField()
    start = serializers.IntegerField()
    end = serializers.IntegerField()
    strand = serializers.CharField(allow_null=True)
    coordinates = serializers.CharField()


class TranscriptVersionDetailSerializer(serializers.ModelSerializer):
    """ A transcript's build specific location. Coordinates come out of the cdot data blob via
        TranscriptVersion.get_coordinates(), which is null when an older import didn't store them """
    accession = serializers.CharField(read_only=True)
    genome_build = serializers.CharField(source="genome_build_id")
    contig = GeneContigSerializer()
    # DRF calls get_coordinates() once and yields null when it returns None (older imports lack coordinates)
    coordinates = GenomicCoordinatesSerializer(source="get_coordinates", read_only=True)

    class Meta:
        model = TranscriptVersion
        fields = ('accession', 'version', 'genome_build', 'biotype', 'contig', 'coordinates')


class GeneVersionDetailSerializer(serializers.ModelSerializer):
    """ A gene's details for one genome build, with the transcripts that place it on the genome """
    accession = serializers.CharField(read_only=True)
    gene_symbol = serializers.CharField(source="gene_symbol_id", allow_null=True)
    genome_build = serializers.CharField(source="genome_build_id")
    contigs = serializers.SerializerMethodField()
    transcript_versions = serializers.SerializerMethodField()

    class Meta:
        model = GeneVersion
        fields = ('accession', 'version', 'gene_symbol', 'genome_build', 'hgnc_identifier',
                  'biotype', 'description', 'contigs', 'transcript_versions')

    @staticmethod
    def _transcript_versions(obj: GeneVersion):
        return obj.transcriptversion_set.order_by("transcript_id", "version").select_related("contig")

    def get_contigs(self, obj: GeneVersion):
        """ Contigs this gene's transcripts are on - usually one, but a gene can also map to alt scaffolds """
        contigs = {tv.contig for tv in self._transcript_versions(obj)}
        return GeneContigSerializer(sorted(contigs, key=lambda c: c.pk), many=True).data

    def get_transcript_versions(self, obj: GeneVersion):
        return TranscriptVersionDetailSerializer(self._transcript_versions(obj), many=True).data


class GeneDetailSerializer(serializers.ModelSerializer):
    """ A Gene and its per genome build versions. Pass 'genome_build' in the serializer context to
        restrict the versions to one build """
    annotation_consortium = serializers.CharField(source="get_annotation_consortium_display")
    versions = serializers.SerializerMethodField()

    class Meta:
        model = Gene
        fields = ('identifier', 'prefixed_identifier', 'annotation_consortium', 'summary', 'versions')

    def get_versions(self, obj: Gene):
        qs = obj.geneversion_set.all()
        if genome_build := self.context.get("genome_build"):
            qs = qs.filter(genome_build=genome_build)
        qs = qs.order_by("genome_build_id", "version").select_related("genome_build")
        return GeneVersionDetailSerializer(qs, many=True).data


class GeneSymbolDetailSerializer(serializers.ModelSerializer):
    """ Everything we know about a gene symbol - its aliases and the genes each consortium assigned to it """
    aliases = serializers.SerializerMethodField()
    genes = serializers.SerializerMethodField()

    class Meta:
        model = GeneSymbol
        fields = ('symbol', 'aliases', 'genes')

    @staticmethod
    def get_aliases(obj: GeneSymbol) -> list[str]:
        """ Other symbols that resolve to the same genes """
        return [s for s in obj.alias_meta.alias_symbol_strs if s != obj.symbol]

    def get_genes(self, obj: GeneSymbol):
        genes = sorted(obj.alias_meta.genes, key=lambda g: g.identifier)
        return GeneDetailSerializer(genes, many=True, context=self.context).data


class GeneInfoSerializer(serializers.ModelSerializer):

    class Meta:
        model = GeneInfo
        fields = ('name', 'description', 'icon_css_class')


class SampleGeneListSerializer(serializers.ModelSerializer):
    active = serializers.SerializerMethodField(read_only=True)
    gene_list = GeneListSerializer()

    class Meta:
        model = SampleGeneList
        fields = ('pk', 'visible', 'gene_list', 'active')

    def get_active(self, obj) -> bool:
        try:
            return obj.sample.activesamplegenelist.sample_gene_list == obj
        except ActiveSampleGeneList.DoesNotExist:
            return False


class GeneCoverageCanonicalTranscriptSerializer(serializers.ModelSerializer):
    transcript_version = TranscriptVersionSerializer()

    class Meta:
        model = GeneCoverageCanonicalTranscript
        exclude = ("gene_coverage_collection", )
        #fields = "__all__"


class GeneCoverageCollectionSerializer(serializers.ModelSerializer):
    genome_build = GenomeBuildSerializer()
    genecoveragecanonicaltranscript_set = GeneCoverageCanonicalTranscriptSerializer(many=True)

    class Meta:
        model = GeneCoverageCollection
        # TODO: Check if "__all__" also does related automatically?
        fields = ["path", "data_state", "genome_build", "genecoveragecanonicaltranscript_set"]

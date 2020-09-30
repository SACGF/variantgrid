from typing import Dict, Any

from django.contrib.auth.models import User
from rest_framework import serializers

from snpdb.models import Variant, Locus, Trio
from snpdb.models.models_clingen_allele import ClinGenAllele
from snpdb.models.models_genome import GenomeBuild, Contig
from snpdb.models.models_variant import Allele, VariantAllele
from snpdb.models.models_vcf import Project
from snpdb.variant_links import variant_link_info


class UserSerializer(serializers.ModelSerializer):

    class Meta:
        model = User
        fields = ('username', 'first_name', 'last_name', 'email', 'is_staff')


class GenomeBuildSerializer(serializers.ModelSerializer):

    class Meta:
        model = GenomeBuild
        fields = ('name', 'alias')


class ContigSerializer(serializers.ModelSerializer):

    class Meta:
        model = Contig
        fields = ('name', 'refseq_accession')


class LocusSerializer(serializers.ModelSerializer):
    contig = ContigSerializer()

    class Meta:
        model = Locus
        fields = '__all__'


class VariantSerializer(serializers.ModelSerializer):
    locus = LocusSerializer()
    alt = serializers.CharField()

    class Meta:
        model = Variant
        fields = '__all__'


class ClinGenAlleleSerializer(serializers.ModelSerializer):
    id = serializers.SerializerMethodField()

    def get_id(self, obj):
        return str(obj)

    class Meta:
        model = ClinGenAllele
        fields = ('id', 'created', 'modified', 'api_response', 'human_url')


class AlleleSerializer(serializers.ModelSerializer):
    clingen_allele = ClinGenAlleleSerializer()
    __str__ = serializers.SerializerMethodField()

    class Meta:
        model = Allele
        fields = ('id', 'clingen_allele', "__str__")

    def get___str__(self, obj):
        return str(obj)


class VariantAlleleSerializer(serializers.ModelSerializer):
    variant = VariantSerializer()
    genome_build = GenomeBuildSerializer()
    allele = AlleleSerializer()
    origin = serializers.SerializerMethodField()
    conversion_tool = serializers.SerializerMethodField()

    class Meta:
        model = VariantAllele
        fields = '__all__'

    def get_origin(self, obj):
        return obj.get_origin_display()

    def get_conversion_tool(self, obj):
        return obj.get_conversion_tool_display()

    @staticmethod
    def data_with_link_data(variant_allele: 'VariantAllele') -> Dict[str, Any]:
        variant_data = VariantAlleleSerializer(variant_allele).data
        link_data = variant_link_info(variant_allele.variant)
        data = dict(variant_data)
        data['link_data'] = link_data
        return data


class TrioSerializer(serializers.ModelSerializer):
    mother = serializers.SerializerMethodField()
    father = serializers.SerializerMethodField()
    proband = serializers.SerializerMethodField()

    class Meta:
        model = Trio
        fields = '__all__'

    def get_mother(self, obj):
        return obj.mother.name

    def get_father(self, obj):
        return obj.father.name

    def get_proband(self, obj):
        return obj.proband.name


class ProjectSerializer(serializers.ModelSerializer):

    class Meta:
        model = Project
        fields = ('name', 'description')

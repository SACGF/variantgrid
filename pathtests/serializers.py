from rest_framework import serializers

from genes.serializers import GeneSerializer, GeneListSerializer
from pathtests.models import PathologyTestVersion, PathologyTest, \
    PathologyTestGeneModificationRequest
from snpdb.serializers import UserSerializer


class PathologyTestSerializer(serializers.ModelSerializer):
    curator = UserSerializer()

    class Meta:
        model = PathologyTest
        fields = '__all__'


class PathologyTestGeneModificationRequestSerializer(serializers.ModelSerializer):
    gene = GeneSerializer()
    user = UserSerializer()

    class Meta:
        model = PathologyTestGeneModificationRequest
        fields = '__all__'


class PathologyTestVersionSerializer(serializers.ModelSerializer):
    pathology_test = PathologyTestSerializer()
    gene_list = GeneListSerializer()
    is_active_test = serializers.SerializerMethodField()
    pathologytestgenemodificationrequest_set = PathologyTestGeneModificationRequestSerializer(many=True, read_only=True)
    __str__ = serializers.SerializerMethodField()

    class Meta:
        model = PathologyTestVersion
        fields = '__all__'

    def get_is_active_test(self, obj):
        active_test = obj.pathology_test.get_active_test_version()
        return obj == active_test

    def get___str__(self, obj):
        return str(obj)

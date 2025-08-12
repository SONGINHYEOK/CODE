from rest_framework import serializers
from .models import Project, Target, Batch

class TargetSerializer(serializers.ModelSerializer):
    class Meta:
        model = Target
        fields = '__all__'

class BatchSerializer(serializers.ModelSerializer):
    target_name = serializers.CharField(source='target.name', read_only=True)
    
    class Meta:
        model = Batch
        fields = '__all__'

class ProjectSerializer(serializers.ModelSerializer):
    leader_name = serializers.CharField(source='leader.username', read_only=True)
    batch_count = serializers.IntegerField(source='batches.count', read_only=True)
    
    class Meta:
        model = Project
        fields = '__all__'
        read_only_fields = ['id', 'created_at', 'updated_at']

class ProjectDetailSerializer(ProjectSerializer):
    batches = BatchSerializer(many=True, read_only=True)
    
    class Meta(ProjectSerializer.Meta):
        fields = '__all__'

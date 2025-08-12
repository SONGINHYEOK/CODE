from rest_framework import viewsets, permissions, status
from rest_framework.decorators import action
from rest_framework.response import Response
from .models import ScreeningJob, ScreeningResult
from .serializers import ScreeningJobSerializer, ScreeningResultSerializer
from .tasks import run_unidock_screening

class ScreeningJobViewSet(viewsets.ModelViewSet):
    queryset = ScreeningJob.objects.all()
    serializer_class = ScreeningJobSerializer
    permission_classes = [permissions.IsAuthenticated]
    
    @action(detail=True, methods=['post'])
    def start(self, request, pk=None):
        job = self.get_object()
        
        if job.status != 'queued':
            return Response(
                {'error': 'Job already started'},
                status=status.HTTP_400_BAD_REQUEST
            )
        
        # Start Celery task
        task = run_unidock_screening.delay(str(job.id))
        job.celery_task_id = task.id
        job.save()
        
        return Response({'task_id': task.id, 'status': 'started'})
    
    @action(detail=True, methods=['get'])
    def results(self, request, pk=None):
        job = self.get_object()
        results = job.results.all()
        serializer = ScreeningResultSerializer(results, many=True)
        return Response(serializer.data)
# omnidrug/urls.py 수정
from django.contrib import admin
from django.urls import path, include
from django.http import JsonResponse
from rest_framework.routers import DefaultRouter
from drf_spectacular.views import SpectacularAPIView, SpectacularSwaggerView
from apps.core.views import CustomTokenObtainPairView, RegisterView, ProfileView
from apps.projects.views import ProjectViewSet, TargetViewSet, BatchViewSet

def home(request):
    return JsonResponse({
        'message': 'Welcome to OMNIDRUG-AI Platform',
        'version': '1.0.0',
        'api_docs': 'http://127.0.0.1:8000/api/docs/',
        'admin': 'http://127.0.0.1:8000/admin/',
        'endpoints': {
            'projects': '/api/projects/',
            'targets': '/api/targets/',
            'batches': '/api/batches/',
            'auth': {
                'login': '/api/auth/login/',
                'register': '/api/auth/register/',
                'profile': '/api/auth/profile/'
            }
        }
    })

router = DefaultRouter()
router.register(r'projects', ProjectViewSet)
router.register(r'targets', TargetViewSet)
router.register(r'batches', BatchViewSet)

urlpatterns = [
    path('', home, name='home'),  # 메인 페이지 추가
    path('admin/', admin.site.urls),
    path('api/', include(router.urls)),
    
    # Authentication
    path('api/auth/login/', CustomTokenObtainPairView.as_view(), name='token_obtain_pair'),
    path('api/auth/register/', RegisterView.as_view(), name='register'),
    path('api/auth/profile/', ProfileView.as_view(), name='profile'),
    
    # API Documentation
    path('api/schema/', SpectacularAPIView.as_view(), name='schema'),
    path('api/docs/', SpectacularSwaggerView.as_view(url_name='schema'), name='swagger-ui'),
]
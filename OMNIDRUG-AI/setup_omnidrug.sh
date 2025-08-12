#!/bin/bash

# OMNIDRUG-AI Platform Auto Setup Script
# This script creates the entire project structure and files

echo "🚀 Starting OMNIDRUG-AI Platform Setup..."

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if we're in the right directory
if [ ! -d "backend" ]; then
    echo -e "${YELLOW}Creating backend directory...${NC}"
    mkdir -p backend
fi

cd backend

# Create Django project
echo -e "${GREEN}Creating Django project...${NC}"
django-admin startproject omnidrug .

# Create apps directory structure
echo -e "${GREEN}Creating app directories...${NC}"
mkdir -p apps/{core,projects,molecules,screening}

# Create __init__.py files for apps
touch apps/__init__.py
touch apps/core/__init__.py
touch apps/projects/__init__.py
touch apps/molecules/__init__.py
touch apps/screening/__init__.py

# Create manage.py if it doesn't exist (Django should create this)
if [ ! -f "manage.py" ]; then
    cat > manage.py << 'EOF'
#!/usr/bin/env python
"""Django's command-line utility for administrative tasks."""
import os
import sys

if __name__ == '__main__':
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'omnidrug.settings')
    try:
        from django.core.management import execute_from_command_line
    except ImportError as exc:
        raise ImportError(
            "Couldn't import Django. Are you sure it's installed and "
            "available on your PYTHONPATH environment variable? Did you "
            "forget to activate a virtual environment?"
        ) from exc
    execute_from_command_line(sys.argv)
EOF
    chmod +x manage.py
fi

# Create requirements.txt
echo -e "${GREEN}Creating requirements.txt...${NC}"
cat > requirements.txt << 'EOF'
# Django
Django==5.0.1
djangorestframework==3.14.0
django-cors-headers==4.3.1
django-environ==0.11.2
django-filter==23.5
drf-spectacular==0.27.0

# Database
psycopg2-binary==2.9.9
django-redis==5.4.0
redis==5.0.1

# Celery
celery==5.3.4
flower==2.0.1

# Scientific
rdkit==2023.9.4
numpy==1.26.3
pandas==2.1.4
scipy==1.11.4
biopython==1.82

# Authentication
djangorestframework-simplejwt==5.3.1

# Utilities
python-dotenv==1.0.0
gunicorn==21.2.0
Pillow==10.2.0
EOF

# Create settings.py
echo -e "${GREEN}Creating settings.py...${NC}"
cat > omnidrug/settings.py << 'EOF'
import os
from pathlib import Path
from datetime import timedelta
import environ

# Environment variables
env = environ.Env()
environ.Env.read_env()

BASE_DIR = Path(__file__).resolve().parent.parent

SECRET_KEY = env('SECRET_KEY', default='django-insecure-change-this-in-production-!@#$%^&*()')
DEBUG = env.bool('DEBUG', default=True)
ALLOWED_HOSTS = env.list('ALLOWED_HOSTS', default=['localhost', '127.0.0.1'])

# Application definition
INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    
    # Third party
    'rest_framework',
    'corsheaders',
    'django_filters',
    'rest_framework_simplejwt',
    'drf_spectacular',
    
    # Our apps
    'apps.core',
    'apps.projects',
    'apps.molecules',
    'apps.screening',
]

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'corsheaders.middleware.CorsMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

ROOT_URLCONF = 'omnidrug.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

WSGI_APPLICATION = 'omnidrug.wsgi.application'

# Database
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': BASE_DIR / 'db.sqlite3',
    }
}

# Password validation
AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]

# Internationalization
LANGUAGE_CODE = 'en-us'
TIME_ZONE = 'UTC'
USE_I18N = True
USE_TZ = True

# Static files
STATIC_URL = 'static/'
STATIC_ROOT = os.path.join(BASE_DIR, 'staticfiles')
MEDIA_URL = '/media/'
MEDIA_ROOT = os.path.join(BASE_DIR, 'media')

# Default primary key field type
DEFAULT_AUTO_FIELD = 'django.db.models.BigAutoField'

# REST Framework
REST_FRAMEWORK = {
    'DEFAULT_AUTHENTICATION_CLASSES': [
        'rest_framework_simplejwt.authentication.JWTAuthentication',
    ],
    'DEFAULT_PERMISSION_CLASSES': [
        'rest_framework.permissions.IsAuthenticated',
    ],
    'DEFAULT_PAGINATION_CLASS': 'rest_framework.pagination.PageNumberPagination',
    'PAGE_SIZE': 20,
    'DEFAULT_SCHEMA_CLASS': 'drf_spectacular.openapi.AutoSchema',
}

# JWT Settings
SIMPLE_JWT = {
    'ACCESS_TOKEN_LIFETIME': timedelta(minutes=60),
    'REFRESH_TOKEN_LIFETIME': timedelta(days=7),
    'ROTATE_REFRESH_TOKENS': True,
}

# CORS
CORS_ALLOWED_ORIGINS = [
    "http://localhost:8080",
    "http://127.0.0.1:8080",
]

# Custom User Model
AUTH_USER_MODEL = 'core.User'

# Celery Configuration
CELERY_BROKER_URL = env('REDIS_URL', default='redis://localhost:6379/0')
CELERY_RESULT_BACKEND = env('REDIS_URL', default='redis://localhost:6379/0')
CELERY_ACCEPT_CONTENT = ['json']
CELERY_TASK_SERIALIZER = 'json'
CELERY_RESULT_SERIALIZER = 'json'
EOF

# Create core app models.py
echo -e "${GREEN}Creating core app...${NC}"
cat > apps/core/models.py << 'EOF'
from django.contrib.auth.models import AbstractUser
from django.db import models
import uuid

class User(AbstractUser):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    organization = models.CharField(max_length=255, blank=True)
    department = models.CharField(max_length=255, blank=True)
    
    class Meta:
        db_table = 'users'

class TimeStampedModel(models.Model):
    """Abstract base model with created/updated timestamps"""
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    
    class Meta:
        abstract = True
EOF

# Create core app admin.py
cat > apps/core/admin.py << 'EOF'
from django.contrib import admin
from django.contrib.auth.admin import UserAdmin
from .models import User

@admin.register(User)
class CustomUserAdmin(UserAdmin):
    list_display = ['username', 'email', 'organization', 'department', 'is_staff']
    fieldsets = UserAdmin.fieldsets + (
        ('Additional Info', {'fields': ('organization', 'department')}),
    )
EOF

# Create core app serializers.py
cat > apps/core/serializers.py << 'EOF'
from rest_framework import serializers
from django.contrib.auth import get_user_model
from rest_framework_simplejwt.serializers import TokenObtainPairSerializer

User = get_user_model()

class UserSerializer(serializers.ModelSerializer):
    class Meta:
        model = User
        fields = ['id', 'username', 'email', 'first_name', 'last_name', 
                  'organization', 'department']
        read_only_fields = ['id']

class CustomTokenObtainPairSerializer(TokenObtainPairSerializer):
    def validate(self, attrs):
        data = super().validate(attrs)
        data['user'] = UserSerializer(self.user).data
        return data
EOF

# Create core app views.py
cat > apps/core/views.py << 'EOF'
from rest_framework import generics, permissions, status
from rest_framework.response import Response
from rest_framework.views import APIView
from rest_framework_simplejwt.views import TokenObtainPairView
from django.contrib.auth import get_user_model
from .serializers import UserSerializer, CustomTokenObtainPairSerializer

User = get_user_model()

class CustomTokenObtainPairView(TokenObtainPairView):
    serializer_class = CustomTokenObtainPairSerializer

class RegisterView(generics.CreateAPIView):
    queryset = User.objects.all()
    serializer_class = UserSerializer
    permission_classes = [permissions.AllowAny]
    
    def create(self, request, *args, **kwargs):
        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        
        user = User.objects.create_user(
            username=serializer.validated_data['username'],
            email=serializer.validated_data['email'],
            password=request.data.get('password', 'defaultpass123'),
            first_name=serializer.validated_data.get('first_name', ''),
            last_name=serializer.validated_data.get('last_name', ''),
            organization=serializer.validated_data.get('organization', ''),
            department=serializer.validated_data.get('department', '')
        )
        
        return Response(UserSerializer(user).data, status=status.HTTP_201_CREATED)

class ProfileView(APIView):
    permission_classes = [permissions.IsAuthenticated]
    
    def get(self, request):
        serializer = UserSerializer(request.user)
        return Response(serializer.data)
EOF

# Create projects app models.py
echo -e "${GREEN}Creating projects app...${NC}"
cat > apps/projects/models.py << 'EOF'
from django.db import models
from django.contrib.auth import get_user_model
from apps.core.models import TimeStampedModel
import uuid

User = get_user_model()

class Project(TimeStampedModel):
    STATUS_CHOICES = [
        ('planning', 'Planning'),
        ('active', 'Active'),
        ('completed', 'Completed'),
        ('archived', 'Archived'),
    ]
    
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    code = models.CharField(max_length=50, unique=True)
    name = models.CharField(max_length=255)
    description = models.TextField(blank=True)
    leader = models.ForeignKey(User, on_delete=models.PROTECT, related_name='led_projects')
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='planning')
    
    therapeutic_area = models.CharField(max_length=100, blank=True)
    target_indication = models.CharField(max_length=255, blank=True)
    
    class Meta:
        db_table = 'projects'
        ordering = ['-created_at']
    
    def __str__(self):
        return f"{self.code}: {self.name}"

class Target(TimeStampedModel):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=255)
    uniprot_id = models.CharField(max_length=20, unique=True, null=True, blank=True)
    pdb_id = models.CharField(max_length=10, blank=True)
    sequence = models.TextField(blank=True)
    
    class Meta:
        db_table = 'targets'
    
    def __str__(self):
        return self.name

class Batch(TimeStampedModel):
    STATUS_CHOICES = [
        ('pending', 'Pending'),
        ('running', 'Running'),
        ('completed', 'Completed'),
        ('failed', 'Failed'),
    ]
    
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    project = models.ForeignKey(Project, on_delete=models.CASCADE, related_name='batches')
    name = models.CharField(max_length=255)
    target = models.ForeignKey(Target, on_delete=models.PROTECT)
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='pending')
    
    # UniDock parameters
    center_x = models.FloatField(null=True, blank=True)
    center_y = models.FloatField(null=True, blank=True)
    center_z = models.FloatField(null=True, blank=True)
    size_x = models.FloatField(default=20.0)
    size_y = models.FloatField(default=20.0)
    size_z = models.FloatField(default=20.0)
    
    class Meta:
        db_table = 'batches'
        unique_together = ['project', 'name']
    
    def __str__(self):
        return f"{self.project.code} - {self.name}"
EOF

# Create projects app admin.py
cat > apps/projects/admin.py << 'EOF'
from django.contrib import admin
from .models import Project, Target, Batch

@admin.register(Project)
class ProjectAdmin(admin.ModelAdmin):
    list_display = ['code', 'name', 'leader', 'status', 'created_at']
    list_filter = ['status', 'created_at']
    search_fields = ['code', 'name', 'description']

@admin.register(Target)
class TargetAdmin(admin.ModelAdmin):
    list_display = ['name', 'uniprot_id', 'pdb_id']
    search_fields = ['name', 'uniprot_id', 'pdb_id']

@admin.register(Batch)
class BatchAdmin(admin.ModelAdmin):
    list_display = ['name', 'project', 'target', 'status', 'created_at']
    list_filter = ['status', 'created_at']
    search_fields = ['name', 'project__name', 'target__name']
EOF

# Create projects app serializers.py
cat > apps/projects/serializers.py << 'EOF'
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
EOF

# Create projects app views.py
cat > apps/projects/views.py << 'EOF'
from rest_framework import viewsets, permissions, status
from rest_framework.decorators import action
from rest_framework.response import Response
from .models import Project, Target, Batch
from .serializers import (ProjectSerializer, ProjectDetailSerializer, 
                          TargetSerializer, BatchSerializer)

class ProjectViewSet(viewsets.ModelViewSet):
    queryset = Project.objects.all()
    serializer_class = ProjectSerializer
    permission_classes = [permissions.IsAuthenticatedOrReadOnly]
    
    def get_serializer_class(self):
        if self.action == 'retrieve':
            return ProjectDetailSerializer
        return ProjectSerializer
    
    def perform_create(self, serializer):
        serializer.save(leader=self.request.user)
    
    @action(detail=True, methods=['post'])
    def create_batch(self, request, pk=None):
        project = self.get_object()
        serializer = BatchSerializer(data=request.data)
        
        if serializer.is_valid():
            serializer.save(project=project)
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

class TargetViewSet(viewsets.ModelViewSet):
    queryset = Target.objects.all()
    serializer_class = TargetSerializer
    permission_classes = [permissions.IsAuthenticatedOrReadOnly]

class BatchViewSet(viewsets.ModelViewSet):
    queryset = Batch.objects.all()
    serializer_class = BatchSerializer
    permission_classes = [permissions.IsAuthenticatedOrReadOnly]
    
    def get_queryset(self):
        queryset = super().get_queryset()
        project_id = self.request.query_params.get('project', None)
        if project_id:
            queryset = queryset.filter(project_id=project_id)
        return queryset
EOF

# Create molecules app models.py
echo -e "${GREEN}Creating molecules app...${NC}"
cat > apps/molecules/models.py << 'EOF'
from django.db import models
from apps.core.models import TimeStampedModel
import uuid

class Molecule(TimeStampedModel):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    smiles = models.TextField()
    canonical_smiles = models.TextField(blank=True)
    inchi = models.TextField(blank=True)
    inchi_key = models.CharField(max_length=27, unique=True, null=True, blank=True)
    
    molecular_weight = models.FloatField(null=True, blank=True)
    logp = models.FloatField(null=True, blank=True)
    hbd = models.IntegerField(null=True, blank=True)
    hba = models.IntegerField(null=True, blank=True)
    
    source = models.CharField(max_length=50, blank=True)
    
    class Meta:
        db_table = 'molecules'
    
    def __str__(self):
        return self.smiles[:50]
EOF

# Create molecules app admin.py
cat > apps/molecules/admin.py << 'EOF'
from django.contrib import admin
from .models import Molecule

@admin.register(Molecule)
class MoleculeAdmin(admin.ModelAdmin):
    list_display = ['id', 'smiles', 'molecular_weight', 'logp', 'source']
    search_fields = ['smiles', 'canonical_smiles', 'inchi_key']
    list_filter = ['source']
EOF

# Create screening app models.py
echo -e "${GREEN}Creating screening app...${NC}"
cat > apps/screening/models.py << 'EOF'
from django.db import models
from apps.core.models import TimeStampedModel
from apps.projects.models import Batch
from apps.molecules.models import Molecule
import uuid

class ScreeningJob(TimeStampedModel):
    STATUS_CHOICES = [
        ('queued', 'Queued'),
        ('running', 'Running'),
        ('completed', 'Completed'),
        ('failed', 'Failed'),
    ]
    
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    batch = models.ForeignKey(Batch, on_delete=models.CASCADE, related_name='screening_jobs')
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='queued')
    
    # Input files
    receptor_file = models.FileField(upload_to='receptors/', null=True, blank=True)
    ligand_file = models.FileField(upload_to='ligands/', null=True, blank=True)
    
    # UniDock parameters
    exhaustiveness = models.IntegerField(default=8)
    num_modes = models.IntegerField(default=9)
    energy_range = models.FloatField(default=3.0)
    
    # Results
    result_file = models.FileField(upload_to='results/', null=True, blank=True)
    
    # Task tracking
    celery_task_id = models.CharField(max_length=255, blank=True)
    started_at = models.DateTimeField(null=True, blank=True)
    completed_at = models.DateTimeField(null=True, blank=True)
    error_message = models.TextField(blank=True)
    
    class Meta:
        db_table = 'screening_jobs'

class ScreeningResult(TimeStampedModel):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    job = models.ForeignKey(ScreeningJob, on_delete=models.CASCADE, related_name='results')
    molecule = models.ForeignKey(Molecule, on_delete=models.CASCADE)
    
    docking_score = models.FloatField()
    rmsd = models.FloatField(null=True, blank=True)
    binding_pose = models.JSONField(default=dict)
    
    is_hit = models.BooleanField(default=False)
    
    class Meta:
        db_table = 'screening_results'
        ordering = ['docking_score']
EOF

# Create screening app admin.py
cat > apps/screening/admin.py << 'EOF'
from django.contrib import admin
from .models import ScreeningJob, ScreeningResult

@admin.register(ScreeningJob)
class ScreeningJobAdmin(admin.ModelAdmin):
    list_display = ['id', 'batch', 'status', 'created_at', 'completed_at']
    list_filter = ['status', 'created_at']
    search_fields = ['batch__name']

@admin.register(ScreeningResult)
class ScreeningResultAdmin(admin.ModelAdmin):
    list_display = ['id', 'job', 'molecule', 'docking_score', 'is_hit']
    list_filter = ['is_hit', 'job']
    search_fields = ['molecule__smiles']
EOF

# Update URLs
echo -e "${GREEN}Updating URLs...${NC}"
cat > omnidrug/urls.py << 'EOF'
from django.contrib import admin
from django.urls import path, include
from rest_framework.routers import DefaultRouter
from drf_spectacular.views import SpectacularAPIView, SpectacularSwaggerView
from apps.core.views import CustomTokenObtainPairView, RegisterView, ProfileView
from apps.projects.views import ProjectViewSet, TargetViewSet, BatchViewSet

router = DefaultRouter()
router.register(r'projects', ProjectViewSet)
router.register(r'targets', TargetViewSet)
router.register(r'batches', BatchViewSet)

urlpatterns = [
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
EOF

# Create Celery configuration
echo -e "${GREEN}Setting up Celery...${NC}"
cat > omnidrug/celery.py << 'EOF'
import os
from celery import Celery

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'omnidrug.settings')

app = Celery('omnidrug')
app.config_from_object('django.conf:settings', namespace='CELERY')
app.autodiscover_tasks()

@app.task(bind=True)
def debug_task(self):
    print(f'Request: {self.request!r}')
EOF

# Update __init__.py for Celery
cat > omnidrug/__init__.py << 'EOF'
from .celery import app as celery_app

__all__ = ('celery_app',)
EOF

# Create .env file
echo -e "${GREEN}Creating .env file...${NC}"
cd ..
cat > .env << 'EOF'
# Django
SECRET_KEY=django-insecure-change-this-in-production-!@#$%^&*()
DEBUG=True
ALLOWED_HOSTS=localhost,127.0.0.1

# Database (using SQLite for initial setup)
# For PostgreSQL, uncomment and configure:
# DB_NAME=omnidrug
# DB_USER=omnidrug_user
# DB_PASSWORD=omnidrug_pass
# DB_HOST=localhost
# DB_PORT=5432

# Redis
REDIS_URL=redis://localhost:6379

# CORS
CORS_ALLOWED_ORIGINS=http://localhost:8080

# UniDock
UNIDOCK_PATH=/opt/unidock
UNIDOCK_GPU_ENABLED=False
EOF

# Create .gitignore
echo -e "${GREEN}Creating .gitignore...${NC}"
cat > .gitignore << 'EOF'
# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
env/
venv/
ENV/
.venv

# Django
*.log
*.pot
*.pyc
local_settings.py
db.sqlite3
media/
staticfiles/

# IDE
.vscode/
.idea/
*.swp
*.swo

# Environment
.env
.env.local

# OS
.DS_Store
Thumbs.db
EOF

echo -e "${GREEN}✅ Setup complete!${NC}"
echo ""
echo "Next steps:"
echo "1. Install dependencies:"
echo "   cd backend && pip install -r requirements.txt"
echo ""
echo "2. Run migrations:"
echo "   python manage.py makemigrations"
echo "   python manage.py migrate"
echo ""
echo "3. Create superuser:"
echo "   python manage.py createsuperuser"
echo ""
echo "4. Run development server:"
echo "   python manage.py runserver"
echo ""
echo "5. Access admin at: http://127.0.0.1:8000/admin/"
echo "6. Access API docs at: http://127.0.0.1:8000/api/docs/"
from django.urls import path
from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('pipeline/', views.pipeline, name='pipeline'),
    path('dataset/', views.dataset, name='dataset'),
    path('qa/', views.qa, name='qa'),
]
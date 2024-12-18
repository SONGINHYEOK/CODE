# drugdevelopment/urls.py
from django.urls import path
from . import views

urlpatterns = [
    path("login/", views.authenticate_user, name="authenticate_user"),
    path("logout/", views.logout_user, name="logout_user"),
    path("", views.get_home, name="home"),
    # path('', TemplateView.as_view(template_name='pipeline.html'), name='home'),  # Add root URL pattern
    path("pipeline/", views.list_projects, name="list_projects"),
    path("pipeline_modal/", views.create_project, name="create_project"),
    path("compound_batch/", views.get_compound_batch, name="get_compound_batch"),
    path("peptide_batch/", views.get_peptide_batch, name="get_peptide_batch"),
    path("create_job/", views.create_job, name="create_job"),
    path("search_operator/", views.search_operators, name="search_operators"),
    path("available_jobs/<int:project_id>/", views.list_available_jobs, name="list_available_jobs"),    
    path("job/start/", views.get_job_start, name="get_job_start"),
    path("job/data/<int:selected_job_id>/", views.get_job_data, name="get_job_data"),
    path("experiment/", views.create_experiment, name="create_experiment"),
    path("docking/results/", views.get_docking_results, name="get_docking_results"),
    path("molecule_detail/", views.get_molecule_detail, name="get_molecule_detail"),
    path("molecules/", views.list_molecules, name="list_molecules"),
    path("datasets/", views.list_datasets, name="list_datasets"),
    path("dataset/create/", views.create_dataset, name="create_dataset"),
    path("dataset/<int:dataset_id>/delete/", views.delete_dataset, name="delete_dataset"),
    path("dataset/detail/", views.get_dataset_detail, name="get_dataset_detail"),
    path("dataset/compound/<int:dataset_compound_id>/delete/", views.delete_dataset_compound, name="delete_dataset_compound"),
    path("assay/start/", views.get_assay_start, name="get_assay_start"),
    path("assay/results/", views.get_assay_results, name="get_assay_results"),
    path("assay/data/", views.create_assay_data, name="create_assay_data"),
]

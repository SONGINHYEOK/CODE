from django.shortcuts import render
from django.http import JsonResponse
from .models import Target, Project

def index(request):
    return render(request, 'drug_discovery/index.html')

def pipeline(request):
    targets = Target.objects.all()  # 모든 타겟 데이터 가져오기
    return render(request, 'drug_discovery/pipeline.html', {
        'targets': targets
    })

def dataset(request):
    return render(request, 'drug_discovery/dataset.html')

def qa(request):
    return render(request, 'drug_discovery/qa.html')

def create_project(request):
    if request.method == "POST":
        try:
            name = request.POST.get('name')
            target_id = request.POST.get('target')
            progress = request.POST.get('progress')
            chem_or_pep = request.POST.get('chem_or_pep')
            
            target = Target.objects.get(id=target_id)
            project = Project.objects.create(
                name=name,
                target=target,
                progress=f"{progress}%",
                chem_or_pep=chem_or_pep
            )
            
            return JsonResponse({
                'success': True,
                'project_id': project.id,
                'target_id': target.id
            })
        except Exception as e:
            return JsonResponse({
                'success': False,
                'error': str(e)
            })

def protein_viewer(request):
    project_id = request.GET.get('project_id')
    target_id = request.GET.get('target_id')
    
    project = Project.objects.get(id=project_id)
    target = Target.objects.get(id=target_id)
    
    context = {
        'project_id': project_id,
        'target_type': target.pdb_id,  # 또는 target과 관련된 PDB ID
        'project': project,
    }
    
    return render(request, 'protein_viewer.html', context)
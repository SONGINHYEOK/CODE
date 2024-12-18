import json
from django.db import transaction
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib.auth import authenticate, login, logout
from .models import Project, Batch, Job, JobType, JobEdge, Target, TargetComment, TargetGene, TargetPathway, TargetPublication, Protein, Peptide, Bioassay, BioassayResult, Compound, Member, ProjectMember, PDB, ProteinCompound, ProteinPeptide, JobCompound, JobPeptide, ResultCompound, Physicochemical, Molecule, Dataset, DatasetCompound,  DatasetCompoundDetail, TargetPDB, JobPDB
from django.http import JsonResponse, HttpResponse, Http404, HttpResponseServerError
from django.views.decorators.csrf import csrf_exempt, csrf_protect
from django.views.decorators.http import require_http_methods
from django.core import serializers
from django.core.serializers.json import DjangoJSONEncoder
from django.core.exceptions import ObjectDoesNotExist
from django.utils import timezone
from django.db import IntegrityError
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem
import os
from django.conf import settings
import logging
import pandas as pd
from django.db.models import Prefetch
from base64 import b64encode
from io import BytesIO
from rdkit.Chem import Draw
logger = logging.getLogger(__name__)

def authenticate_user(request):
    if request.method == "POST":
        username = request.POST.get("username")
        password = request.POST.get("password")
        user = authenticate(request, username=username, password=password)
        if user is not None:
            login(request, user)
            return redirect("home")
        else:
            # Handle invalid login
            pass
    return render(request, "login.html")

def logout_user(request):
    logout(request)
    return redirect("login")

def get_home(request):
    return render(request, "home.html")

def list_projects(request):
    projects = Project.objects.all()
    targets = Target.objects.all()

    context = {
        "projects": projects,
        "targets": targets
    }
    return render(request, "pipeline.html", context)

# def pipeline_view(request):
    # if request.user.is_authenticated:
        # projects = Project.objects.all()  # DB에서 모든 Project 객체를 가져옴
        # return render(request, 'pipeline.html', {'projects': projects})
    # else:
    #     return redirect('login')

@csrf_exempt
def create_project(request):
    if request.method == "POST":
        try:
            name = request.POST.get("name")
            target = request.POST.get("target")
            target_instance, created = Target.objects.get_or_create(target_type=target)
            progress = request.POST.get("progress")
            chem_or_pep = request.POST.get("chem_or_pep")

            project = Project.objects.create(
                name=name,
                target=target_instance,
                progress=progress,
                chem_or_pep=chem_or_pep
            )

            response_data = {
                "id": project.id,
                "name": project.name,
                "target": project.target.target_type,
                "progress": project.progress,
                "chem_or_pep": project.chem_or_pep
            }
            return JsonResponse(response_data)

        except Exception as e:
            print("Error:", str(e))
            return JsonResponse({"error": str(e)}, status=500)
    
    return JsonResponse({"error": "Invalid request method"}, status=400)

def build_network_data(project):
    network_nodes = []
    network_links = []

    network_nodes.append({
        "id": f"project_{project.id}",
        "type": "project",
        "name": project.name
    })

    target = project.target
    network_nodes.append({
        "id": f"target_{target.id}",
        "type": "target",
        "name": target.target_type
    })
    network_links.append({
        "source": f"project_{project.id}",
        "target": f"target_{target.id}"
    })

    pdbs = TargetPDB.objects.filter(target=target).select_related("pdb")
    for pdb in pdbs:
        network_nodes.append({
            "id": f"pdb_{pdb.id}",
            "type": "pdb",
            "name": pdb.original_id
        })
        network_links.append({
            "source": f"target_{target.id}",
            "target": f"pdb_{pdb.id}"
        })
    
    jobs = Job.objects.filter(batch__project=project).select_related("batch", "job_type")
    for job in jobs:
        network_nodes.append({
            "id": f"job_{job.id}",
            "type": "job",
            "name": f"Job {job.id}"
        })

        job_pdbs = JobPDB.objects.filter(job=job)
        for job_pdb in job_pdbs:
            network_links.append({
                "source": f"pdb_{job_pdb.pdb.id}",
                "target": f"job_{job.id}"
            })
    
    job_compounds = JobCompound.objects.filter(job__in=jobs).select_related("job", "compound")
    for jc in job_compounds:
        network_nodes.append({
            "id": f"jobcompound_{jc.id}",
            "type": "job_compound",
            "name": f"Compound {jc.compound.id}"
        })
        network_links.append({
            "source": f"job_{jc.job.id}",
            "target": f"jobcompound_{jc.id}"
        })
    
    return ({
        "network_nodes": network_nodes,
        "network_links": network_links
    })

def get_compound_batch(request):
    project_id = request.GET.get("project_id")
    project = get_object_or_404(Project, id=project_id)

    # Project와 연결된 모든 Batch 객체를 가져옴. 프로젝트 객체도 preload.
    batches = Batch.objects.filter(project=project).select_related("project")

    # Batch 객체들에 속한 모든 Job 객체를 가져옴.
    # Job과 연결된 Batch, JobType, Job의 Operator와 연결된 Member 객체도 preload.
    all_jobs = Job.objects.filter(batch__in=batches).select_related("batch", "job_type", "operator__member")
    job_types = JobType.objects.all()
    
    batch_jobs = []  # 테이블에 표시할 리스트 초기화
    chart_jobs = []  # 노드 그릴 리스트 초기화
    has_target_job = False

    for job in all_jobs:
        # job_type = job.job_type if job.job_type else None
        operator_name = job.operator.member.name if job.operator else "N/A"
        operator_id = job.operator.id if job.operator else None
        start_date = job.start_date if job.start_date else None
        end_date = job.end_date if job.end_date else None

        job_info = {
            "batch_id": job.batch.id,
            "job_id": job.id,
            "job_type": job.job_type.name,
            "operator": operator_name,
            "operator_id": operator_id,
            "status": job.status if not start_date else (10 if start_date and not end_date else job.status),
            "start": start_date.isoformat() if start_date else "N/A",
            "end": end_date.isoformat() if end_date else "N/A",
        }

        if job.job_type.name != "target":
            batch_jobs.append(job_info)
        
        if job.job_type.name == "target" and not has_target_job:
            chart_jobs.append(job_info)
            has_target_job = True
        elif job.job_type.name != "target" and has_target_job:
            chart_jobs.append(job_info)

    batch_jobs.sort(key=lambda x: x["job_id"])

    job_edges = JobEdge.objects.filter(end__in=all_jobs).select_related("start", "end")
    flowchart_links = [{"source": edge.start.id, "target": edge.end.id} for edge in job_edges if edge.start]
    
    network_data = build_network_data(project)

    print("Network Data:", network_data)
    print("Nodes count:", len(network_data['network_nodes']))
    print("Links count:", len(network_data['network_links']))
    
    context = {
        "all_jobs": all_jobs,
        "project_id": project_id,
        "project": project,
        "batch_jobs_json": json.dumps(batch_jobs),
        "chart_jobs_json": json.dumps(chart_jobs),       
        "job_types": job_types,
        "flowchart_links_json": json.dumps(flowchart_links),
        "network_data_json": json.dumps(network_data)
    }

    return render(request, "compound_batch.html", context)

# def peptide_batch_view(request):
#     project_id = request.GET.get("project_id")
#     # project = get_object_or_404(Project, id=project_id)
#     return render(request, "peptide_batch.html", {"project_id": project_id})

def get_peptide_batch(request):
    project_id = request.GET.get("project_id")
    project = get_object_or_404(Project, id=project_id)

    # Project와 연결된 모든 Batch 객체를 가져옴. 프로젝트 객체도 preload.
    batches = Batch.objects.filter(project=project).select_related("project")

    # Batch 객체들에 속한 모든 Job 객체를 가져옴.
    # Job과 연결된 Batch, JobType, Job의 Operator와 연결된 Member 객체도 preload.
    all_jobs = Job.objects.filter(batch__in=batches).select_related("batch", "job_type", "operator__member")
    job_types = JobType.objects.all()
    
    batch_jobs = []  # 테이블에 표시할 리스트 초기화
    chart_jobs = []  # 노드 그릴 리스트 초기화
    has_target_job = False

    for job in all_jobs:
        # job_type = job.job_type if job.job_type else None
        operator_name = job.operator.member.name if job.operator else "N/A"
        operator_id = job.operator.id if job.operator else None
        start_date = job.start_date if job.start_date else None
        end_date = job.end_date if job.end_date else None

        job_info = {
            "batch_id": job.batch.id,
            "job_id": job.id,
            "job_type": job.job_type.name,
            "operator": operator_name,
            "operator_id": operator_id,
            "status": 10 if start_date else "N/A",
            "start": start_date.isoformat() if start_date else "N/A",
            "end": end_date.isoformat() if end_date else "N/A",
        }

        if job.job_type.name != "target":
            batch_jobs.append(job_info)
        
        if (job.job_type.name == "target" and not has_target_job) or start_date != None:
            chart_jobs.append(job_info)
            if job.job_type.name == "target":
                has_target_job = True

    batch_jobs.sort(key=lambda x: x["job_id"])

    job_edges = JobEdge.objects.filter(end__in=all_jobs).select_related("start", "end")
    links = [{"source": edge.start.id, "target": edge.end.id} for edge in job_edges if edge.start]

    context = {
        "all_jobs": all_jobs,
        "project_id": project_id,
        "project": project,
        "batch_jobs_json": json.dumps(batch_jobs),
        "chart_jobs_json": json.dumps(chart_jobs),       
        "job_types": job_types,
        "links_json": json.dumps(links),
    }
    return render(request, "peptide_batch.html", context)


def create_job(request):
    if request.method == "POST":
        project_id = request.POST.get("project_id")
        job_type_id = request.POST.get("jobType")
        operator_name = request.POST.get("operator")
        previous_jobs = request.POST.getlist("previous_job[]")

        project = get_object_or_404(Project, id=project_id)
        job_type = get_object_or_404(JobType, id=job_type_id)
        operator = get_object_or_404(Member, name=operator_name)

        batch, created = Batch.objects.get_or_create(project=project, defaults={"name": f"Batch for {project.name}"})
        project_member, created = ProjectMember.objects.get_or_create(project=project, member=operator)

        if job_type.name == "target":        
            with transaction.atomic():
                # 기존 음수 id 중 가장 작은 값보다 작은 새 음수 id 생성
                existing_negative_ids = Job.objects.filter(id__lt=0).order_by('id').values_list('id', flat=True)  # id가 0보다 작은 모든 Job 객체들을 필터링하고 id를 오름차순으로 정렬해서 id 값만의 리스트 반환
                if existing_negative_ids:
                    new_id = existing_negative_ids[0] - 1  # 가장 작은 음수 id에서 -1을 뺀 값
                else:
                    new_id = -1  # 음수 id가 없으면 -1부터 시작
            
                job = Job.objects.create(id=new_id, batch=batch, job_type=job_type, operator=project_member)
        else:
            job = Job.objects.create(batch=batch, job_type=job_type, operator=project_member)

        for previous_job_id in previous_jobs:
            if previous_job_id == "start":
                JobEdge.objects.create(start=None, end=job)
            else:
                previous_job = get_object_or_404(Job, id=previous_job_id)
                JobEdge.objects.create(start=previous_job, end=job)
        
        all_jobs = Job.objects.filter(batch__project=project).select_related("job_type")
        available_jobs = [{'id': j.id, 'type': j.job_type.name} for j in all_jobs]

        return JsonResponse({
            "batch_id": batch.id,
            "job_id": job.id,
            "job_type": job_type.name,
            "operator": operator.name,
            "previous_jobs": previous_jobs,
            "available_jobs": available_jobs
        })
    
    return JsonResponse({"error": "Invalid request"}, status=400)

def search_operators(request):
    query = request.GET.get("q", "")  # 검색어를 GET 파라미터에서 가져와
    operators = Member.objects.filter(name__icontains=query)  # 대소문자 구분하지 않는 query 포함 문자열 분별
    results = [{"id": operator.id, "name": operator.name} for operator in operators]
    return JsonResponse(results, safe=False)  # data가 dict가 아닌 경우 장고에서는 safe=False 지정해야 함

def list_available_jobs(request, project_id):
    project = get_object_or_404(Project, id=project_id)
    batches = Batch.objects.filter(project=project)
    jobs = Job.objects.filter(batch__in=batches).selected_related('job_type')

    available_jobs = []
    for job in jobs:
        job_info = {
            'id': job.id,
            'type': job.job_type.name
        }
        available_jobs.append(job_info)
    
    return JsonResponse({'jobs': available_jobs})

def get_job_start(request):
    project_id = request.GET.get("project_id")
    # 현재 프로젝트
    current_project = get_object_or_404(Project, id=project_id)
    # 나머지 프로젝트
    # other_projects = Project.objects.exclude(id=project_id)
    job_id = request.GET.get("job_id")
    current_job = get_object_or_404(Job, id=job_id)

    # 현재 프로젝트와 동일한 target_type을 가진 모든 프로젝트를 가져옴
    same_target_projects = Project.objects.filter(target__target_type=current_project.target.target_type)

    # 동일한 target_type을 가진 프로젝트들의 job 중에서 현재 job을 제외한 모든 job을 가져옴
    other_jobs = Job.objects.filter(batch__project__in=same_target_projects).exclude(id=job_id)

    pdbs = PDB.objects.filter(targetpdb__target=current_project.target)
    pdbs_json = json.dumps(list(pdbs.values('id', 'original_id')), cls=DjangoJSONEncoder)

    batch_id = request.GET.get("batch_id")
    operator_id = request.GET.get("operator_id")

    context = {
        "project_id": project_id,
        # "current_project": current_project,
        # "other_projects": other_projects,
        "batch_id": batch_id,
        "job_id": job_id,
        "current_job": current_job,
        "other_jobs": other_jobs,
        "operator_id": operator_id,
        "target_type": current_project.target.target_type,
        "pdbs": pdbs,
        "pdbs_json": pdbs_json
    }
    return render(request, "job_start.html", context)

@csrf_exempt
def get_job_data(request, selected_job_id):  
    selected_job = get_object_or_404(Job, id=selected_job_id)
    job_compounds = JobCompound.objects.filter(job=selected_job)
    compound_data = []
    for job_compound in job_compounds:
        compound_data.append({
            "id": job_compound.compound.id,
            "type": "Chem",
            "data": job_compound.compound.canonical_smiles,
            "mol_3d": job_compound.compound.molecule.mol_3d if job_compound.compound.molecule else None
        })
    
    return JsonResponse(compound_data, safe=False)

    # if project.chem_or_pep == "Chem":
    #     # 현재 project와 연관된 target으로부터 관련된 모든 protein을 가져옴.    
    #     proteins = Protein.objects.filter(target=project.target)

    #     # 가져온 protein과 연관된 모든 compound를 가져옴.
    #     # Compound 모델에서 ProteinCompound 모델의 protein 필드가 proteins 쿼리셋에 포함된 레코드들을 필터링함.
    #     compounds = Compound.objects.filter(proteincompound__protein__in=proteins).distinct()
    #     related_data = [{"id": c.id, "type": "Chem", "data": c.canonical_smiles} for c in compounds]
    # elif project.chem_or_pep == "Pep":
    #     proteins = Protein.objects.filter(target=project.target)
    #     peptides = Peptide.objects.filter(proteinpeptide__protein__in=proteins).distinct()
    #     related_data = [{"id": p.id, "type": "Pep", "data": p.sequence} for p in peptides]
    
    # print(related_data)
    
    # return JsonResponse(related_data, safe=False)
@csrf_exempt
def create_experiment(request):
    if request.method == 'POST':
        logger.info("Received POST request to start_experiment_view")
        
        job_id = request.POST.get('job_id')
        job = get_object_or_404(Job, id=job_id)
        PDB_ID = request.POST.get('PDB_ID')
        binding_sites = request.POST.get('binding_sites')
        selected_compound_ids = json.loads(request.POST.get('selected_compound_ids', '[]'))

        logger.info(f"Received data:")
        logger.info(f"job_id: {job_id}")
        logger.info(f"PDB_ID: {PDB_ID}")
        logger.info(f"binding_sites: {binding_sites}")
        logger.info(f"selected_compound_ids: {selected_compound_ids}")

        try:
            binding_sites = json.loads(binding_sites) if binding_sites else {}
        except json.JSONDecodeError:
            logger.error("Failed to parse binding_sites JSON data")
            return JsonResponse({"error": "Invalid binding_sites JSON data"}, status=400)

        # Directory to store the .sdf files
        save_dir = os.path.join(settings.MEDIA_ROOT, 'sdf_files')
        os.makedirs(save_dir, exist_ok=True)

        # Process molecules
        processed_molecules = []
        for compound_id in selected_compound_ids:
            try:
                molecule = Molecule.objects.get(compound__id=compound_id)
                mol_3d = molecule.mol_3d
                
                if mol_3d:
                    mol = Chem.MolFromMolBlock(mol_3d)
                    if mol:
                        sdf_filename = f"{compound_id}.sdf"
                        sdf_filepath = os.path.join(save_dir, sdf_filename)
                        
                        with open(sdf_filepath, 'w') as sdf_file:
                            sdf_file.write(Chem.MolToMolBlock(mol))
                        
                        processed_molecules.append(compound_id)
                        logger.info(f"Saved SDF file for molecule {compound_id}")
                    else:
                        logger.warning(f"Failed to parse mol_3d for compound {compound_id}")
                else:
                    logger.warning(f"No mol_3d data for compound {compound_id}")
            except Molecule.DoesNotExist:
                logger.error(f"Molecule not found for compound_id: {compound_id}")

        # Update job start date
        job.start_date = timezone.now()
        job.save()

        logger.info(f"Experiment started successfully for job_id: {job_id}")
        logger.info(f"Processed molecules: {processed_molecules}")

        return JsonResponse({
            "success": True, 
            "message": "Experiment started successfully!", 
            "job_id": job_id,
            "processed_molecules": processed_molecules
        })
    else:
        logger.error(f"Invalid request method: {request.method}")
        return JsonResponse({"error": "Invalid request method"}, status=400)
    
# def search_projects(request):
#     # 검색 가능하도록
#     query = request.GET.get('q', '')
#     if query:
#         projects = Project.objects.filter(name__icontains=query)
#     else:
#         projects = Project.objects.all()
    
#     project_list = [{"id": project.id, "name": project.name} for project in projects]
#     return JsonResponse(project_list, safe=False)
def get_docking_results(request):
    project_id = request.GET.get("project_id")
    project = get_object_or_404(Project, id=project_id)
    job_id = request.GET.get("job_id")
    job = get_object_or_404(Job, id=job_id)

    print(f"Job ID: {job_id}")

    job_compounds = JobCompound.objects.filter(job=job)
    print(f"Found {job_compounds.count()} job compounds")
    
    compounds_data = []
    physico_chemical_data = []
    
    for job_compound in job_compounds:
        print(f"\nProcessing job_compound ID: {job_compound.id}")
        
        # ResultCompound 조회 과정 디버깅
        result_compound = ResultCompound.objects.filter(job_compound=job_compound).first()
        print(f"Result compound found: {result_compound is not None}")
        
        if result_compound:
            print(f"Result value: {result_compound.value}")
            
            compound = job_compound.compound
            print(f"Compound found: {compound is not None}")
            
            if compound:
                print(f"Compound ID: {compound.id}")
                print(f"SMILES: {compound.canonical_smiles}")
                
                compounds_data.append({
                    "job_compound_id": job_compound.id,
                    "result_compound_id": result_compound.id,
                    "smiles": compound.canonical_smiles,
                    "score": result_compound.value,
                    "mol_3d": None  # molecule 데이터는 일단 None으로
                })

                if hasattr(compound, 'physicochemical'):
                    print(f"Has physicochemical data: {hasattr(compound, 'physicochemical')}")
                    physico = compound.physicochemical
                    physico_chemical_data.append({
                        "result_compound_id": result_compound.id,
                        "compound_id": compound.id,
                        "score": result_compound.value,
                        "logP": physico.logP,
                        "molecular_weight": physico.molecular_weight,
                        "hydrogen_bond_acceptors": physico.hydrogen_bond_acceptors,
                        "hydrogen_bond_donors": physico.hydrogen_bond_donors,
                        "lipinski": physico.lipinski,
                        "qed": physico.qed,
                        "stereo_centers": physico.stereo_centers,
                        "tpsa": physico.tpsa
                    })
    
    print(f"\nProcessed {len(compounds_data)} compounds")
    print(f"Processed {len(physico_chemical_data)} physicochemical data")

    context = {
        "compounds_json": json.dumps(compounds_data),
        "physico_json": json.dumps(physico_chemical_data)
    }
   
    return render(request, "docking_result.html", context)


# def molecule_list_view(request):    
#     return render(request, "molecule_list.html")

# def molecule_list_view(request):
#     if request.method == "POST":
#         try:
#             selected_items = request.POST.getlist("selected_items")
#             # selected_binding_sites = request.POST.getlist("selected_binding_sites")
#             # job_id = request.POST.get("job_id")
#             # job = get_object_or_404(Job, id=job_id)
#             # target_type = request.POST.get("target_type")

#             if not selected_items:
#                 return JsonResponse({"error": "No items selected"}, status=400)
            
#             for item in selected_items:
#                 print(item)

#         except Exception as e:
#             return JsonResponse({"error": str(e)}, status=500)
#     return JsonResponse({"error": "Invalid request"}, status=400)

def get_molecule_image(smiles, size=(300, 300)):  # 기본 크기를 더 크게 설정
    """SMILES 문자열로부터 base64로 인코딩된 PNG 이미지 생성"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        
        # 더 선명한 이미지를 위해 크기 증가
        fig = Draw.MolToImage(mol, size=size)
        img_buffer = BytesIO()
        fig.save(img_buffer, format='PNG', quality=95)  # 품질 향상
        
        encoded_img = b64encode(img_buffer.getvalue()).decode('utf-8')
        return f"data:image/png;base64,{encoded_img}"
    except Exception as e:
        print(f"Error generating molecule image: {e}")
        return ""

def get_molecule_detail(request):
    try:
        job_compound_id = request.GET.get("job_compound_id")
        result_compound_id = request.GET.get("result_compound_id")
        result_compound = get_object_or_404(ResultCompound, id=result_compound_id)
        molecule_detail_compound = result_compound.job_compound.compound
        physico_chemical = molecule_detail_compound.physicochemical
        
        # 분자 구조 이미지 생성
        molecule_img = ""
        if molecule_detail_compound.canonical_smiles:
            molecule_img = get_molecule_image(molecule_detail_compound.canonical_smiles)
        
        related_jobs = Job.objects.filter(
            jobcompound__compound=molecule_detail_compound
        ).select_related('batch__project')

        job_project_info = [{
            'job_id': job.id,
            'job_type': job.job_type.name if job.job_type else 'No Type',
            'project_id': job.batch.project.id if job.batch and job.batch.project else None,
            'project_name': job.batch.project.name if job.batch and job.batch.project else 'No Project',
            'project_chem_or_pep': job.batch.project.chem_or_pep if job.batch and job.batch.project else 'No Project',
            'target': job.batch.project.target.target_type
        } for job in related_jobs]

        context = {
            "compound_id": molecule_detail_compound.id,
            "smiles": molecule_detail_compound.canonical_smiles,
            "molecule_img": molecule_img,  # SVG 대신 PNG 이미지로 변경
            "logP": physico_chemical.logP,
            "molecular_weight": physico_chemical.molecular_weight,
            "hydrogen_bond_acceptors": physico_chemical.hydrogen_bond_acceptors,
            "hydrogen_bond_donors": physico_chemical.hydrogen_bond_donors,
            "lipinski": physico_chemical.lipinski,
            "qed": physico_chemical.qed,
            "stereo_centers": physico_chemical.stereo_centers,
            "tpsa": physico_chemical.tpsa,
            "job_project_info": job_project_info
        }

        return render(request, "molecule_detail.html", context)
    except Exception as e:
        print(f"Error in molecule_detail_view: {e}")
        raise

def list_molecules(request):
    if request.method == 'POST':
        # JSON 형식의 데이터를 파싱
        selected_compounds_json = request.body.decode('utf-8')
        selected_compounds = json.loads(selected_compounds_json).get('selected_compounds', [])
        print(selected_compounds)
    else:
        selected_compounds = []
        
    context = {
        'selected_compounds': selected_compounds
    }
    
    return render(request, 'molecule_list.html', context)

def get_assay_start(request):
    project_id = request.GET.get("project_id")
    job_id = request.GET.get("job_id")
    print(f"Received job_id: {job_id}")
    current_job = get_object_or_404(Job, id=job_id)
    # 이전 job의 화합물 데이터 가져오기
    previous_jobs = Job.objects.filter(
        job_end_set__start=current_job
    ).select_related('batch', 'job_type')
    compounds_data = []
    if previous_jobs:
        previous_job = previous_jobs[0]
        job_compounds = JobCompound.objects.filter(
            job=previous_job
        ).select_related('compound__molecule')
        for job_compound in job_compounds:
            compounds_data.append({
                "id": job_compound.compound.id,
                "smiles": job_compound.compound.canonical_smiles,
                "mol_3d": job_compound.compound.molecule.mol_3d if job_compound.compound.molecule else None
            })
    context = {
        "project_id": project_id,
        "job_id": job_id,
        "current_job": current_job,
        "compounds_data": json.dumps(compounds_data)
    }
    print(f"Context being sent to template: {context}")
    return render(request, "assay_start.html", context)

@csrf_exempt
def create_assay_data(request):
    if request.method == 'POST':
        print(request)
        try:
            data=json.loads(request.body)
            print('upload assay data', data)
            job_id = data.get('job_id')
            print('upload assay data job_id', job_id)
            experiment_data = data.get('data')
            print('upload assay data experiment_data', experiment_data)

            job = get_object_or_404(Job, id=job_id)
            print('upload assay data job', job)
            if 'file' in request.FILES:
                file = request.FILES['file']
                # xlsx 또는 csv 파일 처리
                if file.name.endswith('.xlsx'):
                    df = pd.read_excel(file)
                elif file.name.endswith('.csv'):
                    df = pd.read_csv(file)
                else:
                    return JsonResponse({"error": "Invalid file format"}, status=400)
                try:
                    with transaction.atomic():
                        for _, row in df.iterrows():
                            compound = get_object_or_404(Compound, canonical_smiles=row['smiles'])
                            job_compound, _ = JobCompound.objects.get_or_create(
                                job=job,
                                compound=compound
                            )
                            result = ResultCompound.objects.create(
                                job_compound=job_compound,
                                value=row['value'],
                                value_type=row['value_type']
                            )
                    return JsonResponse({"success": True})
                except Exception as e:
                    return JsonResponse({"error": str(e)}, status=500)
            else:
                # 직접 입력된 데이터 처리
                try:
                    with transaction.atomic():
                        for item in experiment_data:
                            compound = get_object_or_404(Compound, id=item['compound_id'])
                            job_compound, _ = JobCompound.objects.get_or_create(
                                job=job,
                                compound=compound
                            )
                            result = ResultCompound.objects.create(
                                job_compound=job_compound,
                                value=item['value'],
                                value_type=item['value_type']
                            )
                    return JsonResponse({"success": True})
                except Exception as e:
                    return JsonResponse({"error": str(e)}, status=500)
        except json.JSONDecodeError:
            return JsonResponse({"error": "Invalid JSON data"}, status=400)
        except Exception as e:
            return JsonResponse({"error": str(e)}, status=500)    
    return JsonResponse({"error": "Invalid request"}, status=400)

def get_assay_results(request):
    project_id = request.GET.get("project_id")
    job_id = request.GET.get("job_id")
    print(project_id,job_id)
    job = get_object_or_404(Job, id=job_id)
    job_compounds = JobCompound.objects.filter(job=job).select_related('compound__molecule')
    compounds_data = []
    for job_compound in job_compounds:
        result_compound = ResultCompound.objects.filter(job_compound=job_compound).first()
        compounds_data.append({
            "id": job_compound.compound.id,
            "smiles": job_compound.compound.canonical_smiles,
            "value": result_compound.value if result_compound else None,
            "mol_3d": job_compound.compound.molecule.mol_3d if job_compound.compound.molecule else None
        })
    context = {
        "compounds_json": json.dumps(compounds_data)
    }
    return render(request, "assay_result.html", context)

def list_datasets(request):
    dataset = Dataset.objects.all()

    context = {
        "dataset": dataset
    }

    return render(request, "dataset.html", context)

@require_http_methods(["DELETE"])
def delete_dataset(request, dataset_id):
    try:
        dataset = Dataset.objects.get(id=dataset_id)
        dataset.delete()

        return JsonResponse({
            'status': 'success',
            'message': '데이터셋이 성공적으로 삭제되었습니다.'
        })
    
    except ObjectDoesNotExist:
        return JsonResponse({
            'status': 'error',
            'message': '데이터셋을 찾을 수 없습니다.'
        })
    
    except Exception as e:
        return JsonResponse({
            'status': 'error',
            'message': f'삭제 중 오류가 발생했습니다.: {str(e)}'
        }, status=500)

@require_http_methods(["POST"])
@csrf_protect
def create_dataset(request):
    try:
        name = request.POST.get("name")
        if not name:
            return JsonResponse({"error": "Name field is required."}, status=400)

        dataset = Dataset.objects.create(name=name)

        response_data = {
            "id": dataset.id,
            "name": dataset.name
        }
        return JsonResponse(response_data, status=201)
    
    except Exception as e:
        return JsonResponse({"error": str(e)}, status=500)

def get_dataset_detail(request):
    dataset_id = request.GET.get("dataset_id")
    if not dataset_id:
        raise Http404("Dataset ID is required")

    try:
        dataset = Dataset.objects.prefetch_related(
            Prefetch(
                'datasetcompound_set',
                queryset=DatasetCompound.objects.select_related(
                    'compound',
                    'compound__physicochemical'
                ).prefetch_related(
                    'datasetcompounddetail_set__source_result__job_compound'  # 여기서 result와 job_compound 정보를 가져옴
                )
            )
        ).get(id=dataset_id)

        dataset_compounds_data = []
        for dataset_compound in dataset.datasetcompound_set.all():
            compound = dataset_compound.compound
            if not compound:
                continue

            # DatasetCompoundDetail을 통해 result 정보 가져오기
            detail = dataset_compound.datasetcompounddetail_set.first()  # 첫 번째 detail 사용
            result_info = {}
            if detail and detail.source_result:
                result = detail.source_result
                if result.job_compound:  # job_compound가 있는지 확인
                    result_info = {
                        'job_compound_id': result.job_compound.id,
                        'result_compound_id': result.id
                    }

            compound_data = {
                "dataset_compound_id": dataset_compound.id,
                "compound_id": compound.id,
                "smiles": compound.canonical_smiles,
                "job_compound_id": result_info.get('job_compound_id'),  # result_info에서 가져오기
                "result_compound_id": result_info.get('result_compound_id'),
                "logP": compound.physicochemical.logP if hasattr(compound, 'physicochemical') else None,
                "molecular_weight": compound.physicochemical.molecular_weight if hasattr(compound, 'physicochemical') else None,
            }

            dataset_compounds_data.append(compound_data)

        context = {
            "dataset": dataset,
            "dataset_compounds_json": json.dumps(dataset_compounds_data, default=str),
            "compound_count": len(dataset_compounds_data)
        }

        return render(request, "dataset_detail.html", context)

    except Dataset.DoesNotExist:
        raise Http404("Dataset이 없습니다.")
    except Exception as e:
        logger.error(f"Error in dataset_detail_view: {str(e)}", exc_info=True)
        return HttpResponseServerError("내부 서버 에러 500")
    
@require_http_methods(["DELETE"])
def delete_dataset_compound(request, dataset_compound_id):
    try:
        dataset_compound = DatasetCompound.objects.get(id=dataset_compound_id)
        dataset_compound.delete()

        return JsonResponse({
            'status': 'success',
            'message': '데이터가 성공적으로 삭제되었습니다.'
        })
    
    except ObjectDoesNotExist:
        return JsonResponse({
            'status': 'error',
            'message': '데이터를 찾을 수 없습니다.'
        })

    except Exception as e:
        return JsonResponse({
            'status': 'error',
            'message': f'삭제 중 오류가 발생했습니다: {str(e)}'
        }, status=500)

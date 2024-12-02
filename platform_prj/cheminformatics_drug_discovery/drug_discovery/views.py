from django.shortcuts import render

def index(request):
    return render(request, 'drug_discovery/index.html')

def pipeline(request):
    return render(request, 'drug_discovery/pipeline.html')

def dataset(request):
    return render(request, 'drug_discovery/dataset.html')

def qa(request):
    return render(request, 'drug_discovery/qa.html')
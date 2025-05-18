from systems_biology_framework import SystemsBiologyFramework

framework = SystemsBiologyFramework()
framework.load_network('hierarchical', 'results/hierarchical_network.graphml')

# 예시 1: 약물 억제
framework.apply_perturbation("D01454", 'inhibition', 'hierarchical')
framework.visualize_network("hierarchical_perturbed_inhibition_D01454")

# 예시 2: 모듈 분석
modules = framework.identify_modules("hierarchical", algorithm="louvain")
framework.visualize_modules("hierarchical")

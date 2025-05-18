from systems_biology_framework import SystemsBiologyFramework
import networkx as nx
framework = SystemsBiologyFramework()

network_dir = 'data'

framework.load_network('metabolic', f'{network_dir}/metabolic_network.graphml')
framework.load_network('gene_regulatory', f'{network_dir}/gene_regulatory_network.graphml')
framework.load_network('ppi', f'{network_dir}/ppi_network.graphml')
framework.load_network('drug_target', '/root/systembiology_prj/kegg_drug_target_network.graphml')

framework.create_hierarchical_network(['metabolic', 'gene_regulatory', 'ppi', 'drug_target'])
nx.write_graphml(framework.networks['hierarchical'], "results/hierarchical_network.graphml")
framework.visualize_hierarchical_network()

print("계층적 네트워크 시각화 완료.")

from bioservices import KEGG
import networkx as nx

kegg = KEGG()
kegg.organism = "hsa"  # 인간

# 1. KEGG Drug 리스트 가져오기
drug_ids = kegg.list("drug").split('\n')
drug_ids = [line.split('\t')[0] for line in drug_ids if line.startswith('D')]

# 테스트용 상위 100개만
drug_ids = drug_ids[:100]

G = nx.Graph()

for drug_id in drug_ids:
    try:
        # 2. Drug info 가져오기
        entry = kegg.get(drug_id)
        parsed = kegg.parse(entry)
        name = parsed.get("NAME", [drug_id])[0]
        targets = parsed.get("TARGET", [])
        
        # 3. 약물 노드 추가
        G.add_node(drug_id, type='drug', name=name)
        
        # 4. 타겟 연결
        for t in targets:
            if 'HSA:' in t:
                gene_id = t.split()[0].replace("HSA:", "hsa:")
                G.add_node(gene_id, type='gene')
                G.add_edge(drug_id, gene_id, type='targets')
    except Exception as e:
        print(f"{drug_id} 처리 실패: {str(e)}")

# 5. Graph 저장
nx.write_graphml(G, "kegg_drug_target_network.graphml")

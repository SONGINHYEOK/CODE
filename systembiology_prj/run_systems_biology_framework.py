import os
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from systems_biology_framework import SystemsBiologyFramework

# 데이터 디렉토리 생성
def create_test_data_directory():
    data_dir = "data"
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    return data_dir
# 예시 오믹스 데이터 생성
def create_sample_omics_data(data_dir):
    # 공통 샘플 및 유전자/단백질 ID
    samples = [f"Sample_{i}" for i in range(1, 21)]  # 20개 샘플
    genes = [f"Gene_{i}" for i in range(1, 101)]     # 100개 유전자
    proteins = [f"Protein_{i}" for i in range(1, 81)] # 80개 단백질
    metabolites = [f"Metabolite_{i}" for i in range(1, 51)] # 50개 대사물질
    
    # 무작위 데이터 생성을 위한 시드 설정
    np.random.seed(42)
    
    # 1. 전사체 데이터 생성 (유전자 발현 값)
    gene_expression = np.random.lognormal(mean=1.0, sigma=0.5, size=(len(samples), len(genes)))
    gene_df = pd.DataFrame(gene_expression, index=samples, columns=genes)
    gene_df.to_csv(os.path.join(data_dir, "transcriptome.csv"))
    
    # 2. 단백질체 데이터 생성 (단백질 발현 값)
    protein_expression = np.random.lognormal(mean=0.8, sigma=0.7, size=(len(samples), len(proteins)))
    protein_df = pd.DataFrame(protein_expression, index=samples, columns=proteins)
    protein_df.to_csv(os.path.join(data_dir, "proteome.csv"))
    
    # 3. 대사체 데이터 생성 (대사물질 농도)
    metabolite_conc = np.random.lognormal(mean=0.5, sigma=0.6, size=(len(samples), len(metabolites)))
    metabolite_df = pd.DataFrame(metabolite_conc, index=samples, columns=metabolites)
    metabolite_df.to_csv(os.path.join(data_dir, "metabolome.csv"))
    
    # 4. 샘플 정보 (메타데이터)
    conditions = np.random.choice(['Control', 'Treatment'], size=len(samples))
    metadata = pd.DataFrame({
        'Sample': samples,
        'Condition': conditions,
        'Age': np.random.normal(loc=50, scale=10, size=len(samples)).astype(int),
        'Sex': np.random.choice(['M', 'F'], size=len(samples))
    })
    metadata.to_csv(os.path.join(data_dir, "sample_metadata.csv"), index=False)
    
    print(f"오믹스 데이터 파일이 {data_dir} 디렉토리에 생성되었습니다.")
    
    return {
        'transcriptome': gene_df,
        'proteome': protein_df,
        'metabolome': metabolite_df,
        'metadata': metadata
    }
    
# 예시 시계열 데이터 생성
def create_sample_time_series_data(data_dir):
    # 시간 포인트
    time_points = [0, 2, 4, 8, 12, 24, 48]  # 시간 (시간)
    genes = [f"Gene_{i}" for i in range(1, 31)]  # 30개 유전자
    
    # 데이터 생성
    data = []
    np.random.seed(42)
    
    for gene in genes:
        # 기본 발현 패턴 생성 (사인파 또는 임의 패턴)
        if np.random.rand() < 0.5:
            # 사인파 패턴
            period = np.random.uniform(12, 48)  # 주기: 12~48시간
            phase = np.random.uniform(0, 2*np.pi)  # 위상
            amplitude = np.random.uniform(0.5, 2.0)  # 진폭
            base = np.random.uniform(1.0, 3.0)  # 기준값
            
            expression = base + amplitude * np.sin(2*np.pi*np.array(time_points)/period + phase)
        else:
            # 임의 패턴 (시간에 따른 연속적 변화)
            expression = np.random.normal(1.0, 0.3, size=len(time_points))
            # 연속적인 추세 추가
            for i in range(1, len(expression)):
                expression[i] = 0.7 * expression[i] + 0.3 * expression[i-1]
            
            # 일부 유전자는 증가/감소 추세 추가
            if np.random.rand() < 0.3:
                # 증가 추세
                trend = np.linspace(0, 1, len(time_points))
                expression = expression + trend
            elif np.random.rand() < 0.3:
                # 감소 추세
                trend = np.linspace(1, 0, len(time_points))
                expression = expression + trend
        
        # 약간의 노이즈 추가
        expression += np.random.normal(0, 0.1, size=len(time_points))
        
        # 음수 방지
        expression = np.maximum(0, expression)
        
        # 데이터 추가
        for i, t in enumerate(time_points):
            data.append({'Time': t, 'GeneID': gene, 'Expression': expression[i]})
    
    # 데이터프레임 생성 및 저장
    time_series_df = pd.DataFrame(data)
    time_series_df.to_csv(os.path.join(data_dir, "time_series.csv"), index=False)
    
    print(f"시계열 데이터 파일이 {data_dir} 디렉토리에 생성되었습니다.")
    
    return time_series_df

# 예시 네트워크 생성 및 저장
def create_sample_networks(data_dir):
    # 1. 대사 네트워크 생성
    metabolic_network = nx.DiGraph()
    
    # 대사물질 노드 추가
    metabolites = ["Glucose", "G6P", "F6P", "F1,6BP", "GAP", "DHAP", 
                  "1,3BPG", "3PG", "2PG", "PEP", "Pyruvate", "Lactate",
                  "Acetyl-CoA", "Citrate", "Isocitrate", "a-KG", 
                  "Succinyl-CoA", "Succinate", "Fumarate", "Malate", "OAA"]
    
    for m in metabolites:
        metabolic_network.add_node(m, type="metabolite")
    
    # 효소 노드 추가
    enzymes = ["HK", "PGI", "PFK", "ALDO", "TPI", "GAPDH", "PGK", "PGM", 
              "ENO", "PK", "LDH", "PDH", "CS", "ACO", "IDH", "a-KGDH", 
              "SCS", "SDH", "FH", "MDH"]
    
    for e in enzymes:
        metabolic_network.add_node(e, type="enzyme")
    
    # 글리콜리시스 경로 엣지 추가
    glycolysis_edges = [
        ("Glucose", "HK", {"type": "substrate"}), ("HK", "G6P", {"type": "product"}),
        ("G6P", "PGI", {"type": "substrate"}), ("PGI", "F6P", {"type": "product"}),
        ("F6P", "PFK", {"type": "substrate"}), ("PFK", "F1,6BP", {"type": "product"}),
        ("F1,6BP", "ALDO", {"type": "substrate"}), ("ALDO", "GAP", {"type": "product"}),
        ("ALDO", "DHAP", {"type": "product"}), ("DHAP", "TPI", {"type": "substrate"}),
        ("TPI", "GAP", {"type": "product"}), ("GAP", "GAPDH", {"type": "substrate"}),
        ("GAPDH", "1,3BPG", {"type": "product"}), ("1,3BPG", "PGK", {"type": "substrate"}),
        ("PGK", "3PG", {"type": "product"}), ("3PG", "PGM", {"type": "substrate"}),
        ("PGM", "2PG", {"type": "product"}), ("2PG", "ENO", {"type": "substrate"}),
        ("ENO", "PEP", {"type": "product"}), ("PEP", "PK", {"type": "substrate"}),
        ("PK", "Pyruvate", {"type": "product"}), ("Pyruvate", "LDH", {"type": "substrate"}),
        ("LDH", "Lactate", {"type": "product"})
    ]
    
    metabolic_network.add_edges_from(glycolysis_edges)
    
    # TCA 사이클 엣지 추가
    tca_edges = [
        ("Pyruvate", "PDH", {"type": "substrate"}), ("PDH", "Acetyl-CoA", {"type": "product"}),
        ("Acetyl-CoA", "CS", {"type": "substrate"}), ("OAA", "CS", {"type": "substrate"}),
        ("CS", "Citrate", {"type": "product"}), ("Citrate", "ACO", {"type": "substrate"}),
        ("ACO", "Isocitrate", {"type": "product"}), ("Isocitrate", "IDH", {"type": "substrate"}),
        ("IDH", "a-KG", {"type": "product"}), ("a-KG", "a-KGDH", {"type": "substrate"}),
        ("a-KGDH", "Succinyl-CoA", {"type": "product"}), ("Succinyl-CoA", "SCS", {"type": "substrate"}),
        ("SCS", "Succinate", {"type": "product"}), ("Succinate", "SDH", {"type": "substrate"}),
        ("SDH", "Fumarate", {"type": "product"}), ("Fumarate", "FH", {"type": "substrate"}),
        ("FH", "Malate", {"type": "product"}), ("Malate", "MDH", {"type": "substrate"}),
        ("MDH", "OAA", {"type": "product"})
    ]
    
    metabolic_network.add_edges_from(tca_edges)
    
    # GraphML 형식으로 저장
    nx.write_graphml(metabolic_network, os.path.join(data_dir, "metabolic_network.graphml"))
    
    # 2. 유전자 조절 네트워크 생성
    gene_network = nx.DiGraph()
    
    # 유전자 노드 추가
    genes = [f"Gene_{i}" for i in range(1, 31)]
    for g in genes:
        gene_network.add_node(g, type="gene")
    
    # 전사인자 노드 추가
    tfs = [f"TF_{i}" for i in range(1, 11)]
    for tf in tfs:
        gene_network.add_node(tf, type="transcription_factor")
    
    # 무작위 조절 관계 생성
    np.random.seed(42)
    for tf in tfs:
        # 각 전사인자는 3-8개의 유전자를 조절
        target_genes = np.random.choice(genes, size=np.random.randint(3, 9), replace=False)
        
        for gene in target_genes:
            # 70%는 활성화, 30%는 억제
            if np.random.rand() < 0.7:
                gene_network.add_edge(tf, gene, type="activation", effect="activation")
            else:
                gene_network.add_edge(tf, gene, type="repression", effect="repression")
    
    # GraphML 형식으로 저장
    nx.write_graphml(gene_network, os.path.join(data_dir, "gene_regulatory_network.graphml"))
    
    # 3. 단백질-단백질 상호작용 네트워크 생성
    ppi_network = nx.Graph()
    
    # 단백질 노드 추가
    proteins = [f"Protein_{i}" for i in range(1, 51)]
    for p in proteins:
        ppi_network.add_node(p, type="protein")
    
    # 일부 효소에 해당하는 단백질 노드 추가 (네트워크 통합 테스트용)
    for i, enzyme in enumerate(enzymes):
        protein_name = f"Protein_E{i+1}"
        ppi_network.add_node(protein_name, type="protein", enzyme=enzyme)
    
    # 무작위 상호작용 생성
    np.random.seed(42)
    # 각 단백질은 평균 3-6개의 다른 단백질과 상호작용
    for p in ppi_network.nodes():
        num_interactions = np.random.randint(3, 7)
        potential_partners = [partner for partner in ppi_network.nodes() if partner != p and not ppi_network.has_edge(p, partner)]
        
        if potential_partners:
            partners = np.random.choice(potential_partners, size=min(num_interactions, len(potential_partners)), replace=False)
            
            for partner in partners:
                # 상호작용 가중치 (0.3-1.0)
                weight = 0.3 + 0.7 * np.random.rand()
                ppi_network.add_edge(p, partner, weight=weight, type="interaction")
    
    # GraphML 형식으로 저장
    nx.write_graphml(ppi_network, os.path.join(data_dir, "ppi_network.graphml"))
    
    print(f"네트워크 파일이 {data_dir} 디렉토리에 생성되었습니다.")
    
    return {
        'metabolic': metabolic_network,
        'gene_regulatory': gene_network,
        'ppi': ppi_network
    }

# 메인 함수
def main():
    try:
        # 1. 테스트 데이터 디렉토리 생성
        data_dir = create_test_data_directory()
        
        # 2. 예시 데이터 생성
        print("예시 데이터 생성 중...")
        omics_data = create_sample_omics_data(data_dir)
        time_series_data = create_sample_time_series_data(data_dir)
        networks = create_sample_networks(data_dir)
        
        # 3. 시스템 바이올로지 프레임워크 초기화 및 실행
        print("\n시스템 바이올로지 프레임워크 실행 중...")
        framework = SystemsBiologyFramework()
        
        # 4. 오믹스 데이터 로드 및 통합
        print("\n1. 오믹스 데이터 로드 및 통합 중...")
        try:
            transcriptome_data = framework.load_omics_data(
                'transcriptome', 
                os.path.join(data_dir, "transcriptome.csv"), 
                id_col=0  # 인덱스 열을 ID로 사용
            )
            
            proteome_data = framework.load_omics_data(
                'proteome', 
                os.path.join(data_dir, "proteome.csv"), 
                id_col=0  # 인덱스 열을 ID로 사용
            )
            
            framework.normalize_omics_data('transcriptome')
            framework.normalize_omics_data('proteome')
            
            integrated_data = framework.integrate_omics_data(['transcriptome_normalized', 'proteome_normalized'], method='pca')
            framework.visualize_integrated_data('pca')
        except Exception as e:
            print(f"오믹스 데이터 통합 중 오류: {str(e)}")
        
        # 5. 네트워크 로드 및 계층적 모델링
        print("\n2. 네트워크 로드 및 계층적 모델링 중...")
        try:
            metabolic_network = framework.load_network('metabolic', os.path.join(data_dir, "metabolic_network.graphml"), format='graphml')
            framework.visualize_network('metabolic')
        except Exception as e:
            print(f"대사 네트워크 로드 중 오류: {str(e)}")
            metabolic_network = framework.load_network('metabolic')  # 내장 예시 네트워크 사용
        
        try:
            gene_network = framework.load_network('gene_regulatory', os.path.join(data_dir, "gene_regulatory_network.graphml"), format='graphml')
            framework.visualize_network('gene_regulatory')
        except Exception as e:
            print(f"유전자 네트워크 로드 중 오류: {str(e)}")
            gene_network = framework.load_network('gene_regulatory')  # 내장 예시 네트워크 사용
        
        try:
            ppi_network = framework.load_network('ppi', os.path.join(data_dir, "ppi_network.graphml"), format='graphml')
        except Exception as e:
            print(f"PPI 네트워크 로드 중 오류: {str(e)}")
            ppi_network = framework.load_network('ppi')  # 내장 예시 네트워크 사용
        
        try:
            hierarchical_network = framework.create_hierarchical_network()
            framework.visualize_hierarchical_network()
        except Exception as e:
            print(f"계층적 네트워크 생성 중 오류: {str(e)}")
        
        # 6. 네트워크 모듈 분석
        print("\n3. 네트워크 모듈 분석 중...")
        try:
            modules = framework.identify_modules('metabolic', algorithm='louvain')
            framework.visualize_modules('metabolic')
        except Exception as e:
            print(f"대사 네트워크 모듈 분석 중 오류: {str(e)}")
        
        try:
            gene_modules = framework.identify_modules('gene_regulatory', algorithm='louvain')
            framework.visualize_modules('gene_regulatory')
        except Exception as e:
            print(f"유전자 네트워크 모듈 분석 중 오류: {str(e)}")
        try:
            drug_target_network = framework.load_network(
                'drug_target',
                '/root/systembiology_prj/kegg_drug_target_network.graphml',
                format='graphml'
            )
            framework.visualize_network('drug_target')
        except Exception as e:
            print(f"약물-타겟 네트워크 로딩 중 오류: {str(e)}")
        
        # 7. 시간적 네트워크 역학 분석
        print("\n4. 시간적 네트워크 역학 분석 중...")
        try:
            time_series_data = framework.load_time_series_data(
                os.path.join(data_dir, "time_series.csv"), 
                time_col='Time', 
                id_col='GeneID'
            )
            dynamic_networks = framework.create_dynamic_network('gene_regulatory', threshold=0.5)
            
            simulation_results = framework.simulate_dynamic_network(steps=10)
            framework.visualize_dynamic_network()
            
            patterns = framework.analyze_dynamic_patterns()
            framework.visualize_dynamic_patterns()
        except Exception as e:
            print(f"시간적 네트워크 분석 중 오류: {str(e)}")
        
        # 8. 다중 스케일 모델링
        print("\n5. 다중 스케일 모델링 중...")
        try:
            multiscale_model = framework.create_multiscale_model(['molecular', 'cellular'], integration_method='vertical')
            simulation_results = framework.run_multiscale_simulation(steps=10)
            framework.visualize_multiscale_simulation()
        except Exception as e:
            print(f"다중 스케일 모델링 중 오류: {str(e)}")
        
        print("\n프레임워크 실행이 완료되었습니다. 결과는 'results' 디렉토리에서 확인할 수 있습니다.")
    
    except Exception as e:
        print(f"오류 발생: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
    
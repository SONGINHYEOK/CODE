# systems_biology_framework.py
import pandas as pd
import numpy as np
import networkx as nx
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, SpectralClustering
from community import community_louvain
from scipy.integrate import solve_ivp
import umap
import json
import os

class SystemsBiologyFramework:
    """
    시스템 바이올로지 네트워크 통합 분석을 위한 프레임워크
    """
    def __init__(self):
        self.omics_data = {}  # 오믹스 데이터 저장
        self.networks = {}    # 네트워크 객체 저장
        self.modules = {}     # 모듈 분석 결과 저장
        self.dynamics = {}    # 시간적 역학 분석 결과 저장
        self.multiscale = {}  # 다중 스케일 모델 저장
        
        # 데이터 및 결과 저장 디렉토리
        self.data_dir = "data"
        self.results_dir = "results"
        
        # 디렉토리 생성
        for dir_path in [self.data_dir, self.results_dir]:
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)
    
    #----------------------------------------------------
    # 1. 다중 오믹스 데이터 통합
    #----------------------------------------------------
    
    def load_omics_data(self, omics_type, file_path, sample_col=None, id_col=None, transpose=False):
        """
        오믹스 데이터 로드
        
        Parameters:
        -----------
        omics_type : str
            데이터 유형 (예: 'transcriptome', 'proteome', 'metabolome', 'genome')
        file_path : str
            데이터 파일 경로 (CSV, TSV 등)
        sample_col : str, optional
            샘플 ID가 포함된 열 이름 (None이면 인덱스 사용)
        id_col : str or int, optional
            특성 ID가 포함된 열 이름 또는 인덱스 (None이면 첫 번째 열 사용)
        transpose : bool, optional
            데이터를 전치할지 여부 (행이 샘플, 열이 특성이 되도록)
        """
        # 데이터 로드
        if file_path.endswith('.csv'):
            if isinstance(id_col, int):
                df = pd.read_csv(file_path, index_col=id_col)
            else:
                df = pd.read_csv(file_path)
        elif file_path.endswith('.tsv') or file_path.endswith('.txt'):
            if isinstance(id_col, int):
                df = pd.read_csv(file_path, sep='\t', index_col=id_col)
            else:
                df = pd.read_csv(file_path, sep='\t')
        elif file_path.endswith('.xlsx'):
            if isinstance(id_col, int):
                df = pd.read_excel(file_path, index_col=id_col)
            else:
                df = pd.read_excel(file_path)
        else:
            raise ValueError("지원되지 않는 파일 형식입니다.")
        
        # 샘플 열 처리
        if sample_col is not None:
            if sample_col in df.columns:
                df = df.set_index(sample_col)
            else:
                raise ValueError(f"샘플 열 '{sample_col}'을 찾을 수 없습니다.")
        
        # ID 열 처리 (문자열인 경우)
        if id_col is not None and isinstance(id_col, str):
            if id_col in df.columns:
                df = df.set_index(id_col)
            else:
                raise ValueError(f"ID 열 '{id_col}'을 찾을 수 없습니다.")
        
        # 데이터 전치 (필요시)
        if transpose:
            df = df.T
        
        # 숫자 데이터만 포함되어 있는지 확인
        numeric_df = pd.DataFrame()
        for col in df.columns:
            try:
                numeric_df[col] = pd.to_numeric(df[col], errors='coerce')
            except:
                print(f"경고: 열 '{col}'에 숫자가 아닌 값이 있습니다. 이 열은 제외됩니다.")
        
        if numeric_df.empty:
            raise ValueError("숫자 데이터를 포함하는 열이 없습니다.")
        
        # 결측치 정보 출력
        print(f"원본 행 수: {numeric_df.shape[0]}, 열 수: {numeric_df.shape[1]}")
        print(f"결측치 비율: {numeric_df.isna().sum().sum() / (numeric_df.shape[0] * numeric_df.shape[1]):.2%}")
        
        # NaN이 50% 이상인 행/열 제거
        numeric_df = numeric_df.dropna(axis=1, thresh=numeric_df.shape[0]*0.5)
        numeric_df = numeric_df.dropna(axis=0, thresh=numeric_df.shape[1]*0.5)
        
        # 나머지 결측치는 평균값으로 대체
        numeric_df = numeric_df.fillna(numeric_df.mean())
        
        print(f"처리 후 행 수: {numeric_df.shape[0]}, 열 수: {numeric_df.shape[1]}")
        
        # 오믹스 데이터 저장
        self.omics_data[omics_type] = numeric_df
        
        return numeric_df
    
    def normalize_omics_data(self, omics_type, method='zscore'):
        """
        오믹스 데이터 정규화
        
        Parameters:
        -----------
        omics_type : str
            정규화할 오믹스 데이터 유형
        method : str, optional
            정규화 방법 ('zscore', 'minmax', 'log', 'quantile')
        """
        if omics_type not in self.omics_data:
            raise ValueError(f"'{omics_type}' 데이터를 찾을 수 없습니다. 먼저 데이터를 로드하세요.")
        
        df = self.omics_data[omics_type].copy()
        
        # 정규화 방법 적용
        if method == 'zscore':
            # Z-score 정규화 (평균 0, 표준편차 1)
            normalized_df = (df - df.mean()) / df.std()
        
        elif method == 'minmax':
            # Min-Max 정규화 (0-1 범위)
            normalized_df = (df - df.min()) / (df.max() - df.min())
        
        elif method == 'log':
            # 로그 변환 (음수 또는 0 값이 있는 경우 처리)
            min_val = df.min().min()
            if min_val <= 0:
                df = df - min_val + 1  # 모든 값을 양수로 변환
            normalized_df = np.log2(df)
        
        elif method == 'quantile':
            # 분위수 정규화
            from sklearn.preprocessing import QuantileTransformer
            qt = QuantileTransformer(output_distribution='normal')
            normalized_values = qt.fit_transform(df.values)
            normalized_df = pd.DataFrame(normalized_values, index=df.index, columns=df.columns)
        
        else:
            raise ValueError(f"지원되지 않는 정규화 방법: {method}")
        
        # 정규화된 데이터 저장
        self.omics_data[f"{omics_type}_normalized"] = normalized_df
        
        return normalized_df
    
    def integrate_omics_data(self, omics_list, method='concatenate', common_samples=True):
        """
        여러 오믹스 데이터셋 통합
        
        Parameters:
        -----------
        omics_list : list
            통합할 오믹스 데이터 유형 목록
        method : str, optional
            통합 방법 ('concatenate', 'correlation', 'pca', 'umap')
        common_samples : bool, optional
            공통 샘플만 사용할지 여부
        """
        # 오믹스 데이터 존재 확인
        for omics in omics_list:
            if omics not in self.omics_data:
                raise ValueError(f"'{omics}' 데이터를 찾을 수 없습니다.")
        
        # 공통 샘플 찾기
        if common_samples:
            common_idx = set(self.omics_data[omics_list[0]].index)
            for omics in omics_list[1:]:
                common_idx = common_idx.intersection(set(self.omics_data[omics].index))
            
            common_idx = sorted(list(common_idx))
            print(f"공통 샘플 수: {len(common_idx)}")
            
            if len(common_idx) == 0:
                raise ValueError("공통 샘플이 없습니다.")
        
        # 통합 방법 적용
        if method == 'concatenate':
            # 특성 단순 연결
            dfs = []
            for omics in omics_list:
                df = self.omics_data[omics]
                if common_samples:
                    df = df.loc[common_idx]
                dfs.append(df)
            
            # 열 이름 충돌 방지
            for i, df in enumerate(dfs):
                df.columns = [f"{omics_list[i]}_{col}" for col in df.columns]
            
            # 데이터프레임 연결
            concatenated_df = pd.concat(dfs, axis=1)
            self.omics_data['integrated'] = concatenated_df
            return concatenated_df
        
        elif method == 'correlation':
            # 오믹스 데이터 간 상관관계 분석
            correlation_matrices = {}
            
            # 각 오믹스 쌍에 대한 상관관계 계산
            for i, omics1 in enumerate(omics_list):
                for j, omics2 in enumerate(omics_list[i+1:], i+1):
                    df1 = self.omics_data[omics1]
                    df2 = self.omics_data[omics2]
                    
                    if common_samples:
                        df1 = df1.loc[common_idx]
                        df2 = df2.loc[common_idx]
                    
                    # 공통 샘플에 대해서만 상관관계 계산
                    common_samples_ij = sorted(set(df1.index).intersection(set(df2.index)))
                    
                    if len(common_samples_ij) < 2:
                        print(f"경고: {omics1}와 {omics2} 간에 공통 샘플이 부족합니다.")
                        continue
                    
                    df1 = df1.loc[common_samples_ij]
                    df2 = df2.loc[common_samples_ij]
                    
                    # 상관관계 행렬 계산 (특성 x 특성)
                    corr_matrix = pd.DataFrame(
                        np.corrcoef(df1.values.T, df2.values.T),
                        index=list(df1.columns) + list(df2.columns),
                        columns=list(df1.columns) + list(df2.columns)
                    )
                    
                    # 결과 저장
                    pair_name = f"{omics1}_vs_{omics2}"
                    correlation_matrices[pair_name] = corr_matrix
            
            self.omics_data['correlation_matrices'] = correlation_matrices
            return correlation_matrices
        
        elif method == 'pca':
            # 각 오믹스 데이터에 PCA 적용 후 주성분 결합
            pca_results = {}
            combined_pcs = []
            
            for omics in omics_list:
                df = self.omics_data[omics]
                if common_samples:
                    df = df.loc[common_idx]
                
                # PCA 적용
                n_components = min(5, df.shape[1])  # 최대 5개 주성분
                pca = PCA(n_components=n_components)
                pca_result = pca.fit_transform(df)
                
                # 결과 저장
                pca_df = pd.DataFrame(
                    pca_result, 
                    index=df.index,
                    columns=[f"{omics}_PC{i+1}" for i in range(n_components)]
                )
                pca_results[omics] = {
                    'pca_df': pca_df,
                    'explained_variance_ratio': pca.explained_variance_ratio_
                }
                
                combined_pcs.append(pca_df)
            
            # 주성분 결합
            combined_pca_df = pd.concat(combined_pcs, axis=1)
            self.omics_data['pca_integrated'] = combined_pca_df
            self.omics_data['pca_results'] = pca_results
            
            return combined_pca_df
        
        elif method == 'umap':
            # 각 오믹스 데이터에 UMAP 적용 후 결합
            umap_results = {}
            combined_embeddings = []
            
            for omics in omics_list:
                df = self.omics_data[omics]
                if common_samples:
                    df = df.loc[common_idx]
                
                # UMAP 적용
                reducer = umap.UMAP(n_components=2, random_state=42)
                embedding = reducer.fit_transform(df)
                
                # 결과 저장
                umap_df = pd.DataFrame(
                    embedding, 
                    index=df.index,
                    columns=[f"{omics}_UMAP1", f"{omics}_UMAP2"]
                )
                umap_results[omics] = umap_df
                combined_embeddings.append(umap_df)
            
            # UMAP 임베딩 결합
            combined_umap_df = pd.concat(combined_embeddings, axis=1)
            self.omics_data['umap_integrated'] = combined_umap_df
            self.omics_data['umap_results'] = umap_results
            
            return combined_umap_df
        
        else:
            raise ValueError(f"지원되지 않는 통합 방법: {method}")
    
    def visualize_integrated_data(self, method='pca', colored_by=None):
        """
        통합된 오믹스 데이터 시각화
        
        Parameters:
        -----------
        method : str, optional
            시각화 방법 ('pca', 'umap', 'heatmap', 'correlation')
        colored_by : str, optional
            색상 구분 기준 (샘플 그룹 등)
        """
        if method == 'pca':
            if 'pca_integrated' not in self.omics_data:
                raise ValueError("먼저 'pca' 방법으로 데이터를 통합하세요.")
            
            pca_df = self.omics_data['pca_integrated']
            
            # 첫 번째와 두 번째 주성분으로 그림
            plt.figure(figsize=(10, 8))
            if colored_by is not None and colored_by in pca_df.index.names:
                # 그룹별 색상
                for group in pca_df.index.get_level_values(colored_by).unique():
                    group_data = pca_df.xs(group, level=colored_by)
                    plt.scatter(
                        group_data.iloc[:, 0], 
                        group_data.iloc[:, 1],
                        label=group,
                        alpha=0.7
                    )
                plt.legend()
            else:
                # 단일 색상
                plt.scatter(pca_df.iloc[:, 0], pca_df.iloc[:, 1], alpha=0.7)
                
                # 샘플 레이블 추가
                for i, sample in enumerate(pca_df.index):
                    plt.annotate(
                        sample, 
                        (pca_df.iloc[i, 0], pca_df.iloc[i, 1]),
                        fontsize=8
                    )
            
            plt.xlabel(pca_df.columns[0])
            plt.ylabel(pca_df.columns[1])
            plt.title("PCA of Integrated Omics Data")
            plt.grid(True, linestyle='--', alpha=0.7)
            
            # 결과 저장
            plt.savefig(os.path.join(self.results_dir, "pca_integrated.png"), dpi=300, bbox_inches='tight')
            plt.close()
        
        elif method == 'umap':
            if 'umap_integrated' not in self.omics_data:
                raise ValueError("먼저 'umap' 방법으로 데이터를 통합하세요.")
            
            umap_df = self.omics_data['umap_integrated']
            
            # UMAP 시각화 (각 오믹스 데이터별)
            num_omics = len(self.omics_data['umap_results'])
            if num_omics > 1:
                fig, axes = plt.subplots(1, num_omics, figsize=(num_omics*5, 5))
                
                for i, (omics, df) in enumerate(self.omics_data['umap_results'].items()):
                    ax = axes[i]
                    ax.scatter(df.iloc[:, 0], df.iloc[:, 1], alpha=0.7)
                    ax.set_title(f"UMAP of {omics} Data")
                    ax.set_xlabel(df.columns[0])
                    ax.set_ylabel(df.columns[1])
                    ax.grid(True, linestyle='--', alpha=0.7)
                
                plt.tight_layout()
            else:
                # 단일 오믹스인 경우
                omics = list(self.omics_data['umap_results'].keys())[0]
                df = self.omics_data['umap_results'][omics]
                
                plt.figure(figsize=(8, 6))
                plt.scatter(df.iloc[:, 0], df.iloc[:, 1], alpha=0.7)
                plt.title(f"UMAP of {omics} Data")
                plt.xlabel(df.columns[0])
                plt.ylabel(df.columns[1])
                plt.grid(True, linestyle='--', alpha=0.7)
            
            # 결과 저장
            plt.savefig(os.path.join(self.results_dir, "umap_integrated.png"), dpi=300, bbox_inches='tight')
            plt.close()
        
        elif method == 'heatmap':
            if 'integrated' not in self.omics_data:
                raise ValueError("먼저 'concatenate' 방법으로 데이터를 통합하세요.")
            
            integrated_df = self.omics_data['integrated']
            
            # 데이터 크기가 너무 크면 샘플링
            if integrated_df.shape[1] > 100:
                print(f"경고: 특성이 너무 많습니다. 처음 100개만 시각화합니다.")
                integrated_df = integrated_df.iloc[:, :100]
            
            if integrated_df.shape[0] > 50:
                print(f"경고: 샘플이 너무 많습니다. 처음 50개만 시각화합니다.")
                integrated_df = integrated_df.iloc[:50, :]
            
            # 히트맵 생성
            plt.figure(figsize=(12, 10))
            sns.heatmap(
                integrated_df,
                cmap='viridis',
                xticklabels=False,  # 특성 이름 숨김 (너무 많을 경우)
                yticklabels=True,   # 샘플 이름 표시
                cbar_kws={'label': 'Normalized Value'}
            )
            plt.title("Heatmap of Integrated Omics Data")
            plt.tight_layout()
            
            # 결과 저장
            plt.savefig(os.path.join(self.results_dir, "heatmap_integrated.png"), dpi=300, bbox_inches='tight')
            plt.close()
        
        elif method == 'correlation':
            if 'correlation_matrices' not in self.omics_data:
                raise ValueError("먼저 'correlation' 방법으로 데이터를 통합하세요.")
            
            corr_matrices = self.omics_data['correlation_matrices']
            
            # 각 상관관계 행렬 시각화
            for pair_name, corr_matrix in corr_matrices.items():
                # 행렬이 너무 크면 샘플링
                if corr_matrix.shape[0] > 100:
                    print(f"경고: {pair_name} 행렬이 너무 큽니다. 처음 100x100만 시각화합니다.")
                    corr_matrix = corr_matrix.iloc[:100, :100]
                
                plt.figure(figsize=(10, 8))
                sns.heatmap(
                    corr_matrix,
                    cmap='coolwarm',
                    vmin=-1,
                    vmax=1,
                    xticklabels=False,
                    yticklabels=False,
                    cbar_kws={'label': 'Correlation Coefficient'}
                )
                plt.title(f"Correlation Heatmap: {pair_name}")
                plt.tight_layout()
                
                # 결과 저장
                plt.savefig(os.path.join(self.results_dir, f"correlation_{pair_name}.png"), dpi=300, bbox_inches='tight')
                plt.close()
        else:
            raise ValueError(f"지원되지 않는 시각화 방법: {method}")
    
    #----------------------------------------------------
    # 2. 계층적 네트워크 모델링
    #----------------------------------------------------
    
    def load_network(self, network_type, file_path=None, format='graphml'):
        """
        네트워크 데이터 로드
        
        Parameters:
        -----------
        network_type : str
            네트워크 유형 (예: 'metabolic', 'gene_regulatory', 'signaling', 'ppi')
        file_path : str, optional
            네트워크 파일 경로
        format : str, optional
            파일 형식 ('graphml', 'gml', 'edgelist', 'adjacency')
        """
        if file_path is None:
            # 예시 네트워크 생성 (실제 구현시 제거 필요)
            G = nx.DiGraph() if network_type in ['metabolic', 'gene_regulatory', 'signaling'] else nx.Graph()
            
            # 예시 노드 추가
            if network_type == 'metabolic':
                G.add_nodes_from(['Glucose', 'G6P', 'F6P', 'FBP', 'PEP', 'Pyruvate'], type='metabolite')
                G.add_nodes_from(['HK', 'PGI', 'PFK', 'ENO', 'PYK'], type='enzyme')
                
                # 예시 엣지 추가
                G.add_edges_from([
                    ('Glucose', 'HK'), ('HK', 'G6P'),
                    ('G6P', 'PGI'), ('PGI', 'F6P'),
                    ('F6P', 'PFK'), ('PFK', 'FBP'),
                    ('PEP', 'PYK'), ('PYK', 'Pyruvate')
                ])
            
            elif network_type == 'gene_regulatory':
                G.add_nodes_from(['Gene1', 'Gene2', 'Gene3', 'Gene4', 'Gene5'], type='gene')
                G.add_nodes_from(['TF1', 'TF2', 'TF3'], type='transcription_factor')
                
                G.add_edges_from([
                    ('TF1', 'Gene1', {'effect': 'activation'}),
                    ('TF1', 'Gene2', {'effect': 'repression'}),
                    ('TF2', 'Gene3', {'effect': 'activation'}),
                    ('TF3', 'Gene4', {'effect': 'activation'}),
                    ('TF3', 'Gene5', {'effect': 'repression'})
                ])
            
            elif network_type == 'signaling':
                G.add_nodes_from(['Receptor1', 'Receptor2'], type='receptor')
                G.add_nodes_from(['Kinase1', 'Kinase2', 'Kinase3'], type='kinase')
                G.add_nodes_from(['TF1', 'TF2'], type='transcription_factor')
                
                G.add_edges_from([
                    ('Receptor1', 'Kinase1', {'effect': 'phosphorylation'}),
                    ('Kinase1', 'Kinase2', {'effect': 'phosphorylation'}),
                    ('Kinase2', 'TF1', {'effect': 'phosphorylation'}),
                    ('Receptor2', 'Kinase3', {'effect': 'phosphorylation'}),
                    ('Kinase3', 'TF2', {'effect': 'phosphorylation'})
                ])
            
            elif network_type == 'ppi':
                G.add_nodes_from(['Protein1', 'Protein2', 'Protein3', 'Protein4', 'Protein5'], type='protein')
                
                G.add_edges_from([
                    ('Protein1', 'Protein2', {'weight': 0.8}),
                    ('Protein1', 'Protein3', {'weight': 0.6}),
                    ('Protein2', 'Protein4', {'weight': 0.7}),
                    ('Protein3', 'Protein4', {'weight': 0.9}),
                    ('Protein4', 'Protein5', {'weight': 0.5})
                ])
            
            else:
                G = nx.DiGraph()  # 기본 빈 네트워크
            
            self.networks[network_type] = G
            return G
        
        # 파일에서 네트워크 로드
        if format == 'graphml':
            G = nx.read_graphml(file_path)
        elif format == 'gml':
            G = nx.read_gml(file_path)
        elif format == 'edgelist':
            G = nx.read_edgelist(file_path, create_using=nx.DiGraph() if network_type in ['metabolic', 'gene_regulatory', 'signaling'] else nx.Graph())
        elif format == 'adjacency':
            # 인접 행렬에서 로드
            adj_matrix = pd.read_csv(file_path, index_col=0)
            G = nx.from_pandas_adjacency(adj_matrix, create_using=nx.DiGraph() if network_type in ['metabolic', 'gene_regulatory', 'signaling'] else nx.Graph())
        else:
            raise ValueError(f"지원되지 않는 파일 형식: {format}")
        
        self.networks[network_type] = G
        return G
    
    def visualize_network(self, network_type, layout='spring', node_color_by='type', node_size_by=None):
        """
        네트워크 시각화
        
        Parameters:
        -----------
        network_type : str
            시각화할 네트워크 유형
        layout : str, optional
            레이아웃 알고리즘 ('spring', 'circular', 'kamada_kawai', 'spectral')
        node_color_by : str, optional
            노드 색상 기준 속성
        node_size_by : str, optional
            노드 크기 기준 속성
        """
        if network_type not in self.networks:
            raise ValueError(f"'{network_type}' 네트워크를 찾을 수 없습니다.")
        
        G = self.networks[network_type]
        
        # 레이아웃 계산
        if layout == 'spring':
            pos = nx.spring_layout(G, seed=42)
        elif layout == 'circular':
            pos = nx.circular_layout(G)
        elif layout == 'kamada_kawai':
            pos = nx.kamada_kawai_layout(G)
        elif layout == 'spectral':
            pos = nx.spectral_layout(G)
        else:
            pos = nx.spring_layout(G, seed=42)  # 기본값
        
        plt.figure(figsize=(12, 10))
        
        # 노드 색상 설정
        if node_color_by and node_color_by in G.nodes[list(G.nodes)[0]]:
            # 노드 유형별 색상 매핑
            node_types = {data[node_color_by] for _, data in G.nodes(data=True) if node_color_by in data}
            color_map = {t: plt.cm.tab10(i) for i, t in enumerate(node_types)}
            
            node_colors = [color_map.get(G.nodes[n].get(node_color_by, 'unknown'), 'gray') for n in G.nodes]
            
            # 범례 추가
            legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, 
                                         label=node_type, markersize=10) 
                             for node_type, color in color_map.items()]
            plt.legend(handles=legend_elements, loc='upper right')
        else:
            node_colors = 'skyblue'
        
        # 노드 크기 설정
        if node_size_by and node_size_by in G.nodes[list(G.nodes)[0]]:
            node_sizes = [300 * G.nodes[n].get(node_size_by, 1) for n in G.nodes]
        else:
            node_sizes = 300
        
        # 엣지 색상 및 스타일 설정
        # 엣지 색상 및 스타일 설정
        edge_colors = []
        edge_styles = []
        
        for u, v, data in G.edges(data=True):
            # 엣지 효과에 따른 색상 설정
            if 'effect' in data:
                if data['effect'] == 'activation' or data['effect'] == 'phosphorylation':
                    edge_colors.append('green')
                elif data['effect'] == 'repression' or data['effect'] == 'inhibition':
                    edge_colors.append('red')
                else:
                    edge_colors.append('gray')
            elif 'weight' in data:
                # 가중치에 따른 색상 설정 (가중치가 높을수록 더 진한 색)
                weight = data['weight']
                edge_colors.append(plt.cm.Blues(0.5 + 0.5 * weight))
            else:
                edge_colors.append('gray')
            
            # 엣지 스타일 설정
            if 'effect' in data and data['effect'] == 'repression':
                edge_styles.append('dashed')
            else:
                edge_styles.append('solid')
        
        # 그래프 그리기
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.8)
        
        # 엣지 스타일에 따라 분리하여 그리기
        solid_edges = [(u, v) for (u, v, d), style in zip(G.edges(data=True), edge_styles) if style == 'solid']
        dashed_edges = [(u, v) for (u, v, d), style in zip(G.edges(data=True), edge_styles) if style == 'dashed']
        
        solid_colors = [color for color, style in zip(edge_colors, edge_styles) if style == 'solid']
        dashed_colors = [color for color, style in zip(edge_colors, edge_styles) if style == 'dashed']
        
        nx.draw_networkx_edges(G, pos, edgelist=solid_edges, edge_color=solid_colors, width=1.5, alpha=0.7)
        nx.draw_networkx_edges(G, pos, edgelist=dashed_edges, edge_color=dashed_colors, 
                              width=1.5, alpha=0.7, style='dashed')
        
        # 노드 레이블 그리기
        nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold')
        
        plt.title(f"{network_type.capitalize()} Network")
        plt.axis('off')
        plt.tight_layout()
        
        # 결과 저장
        plt.savefig(os.path.join(self.results_dir, f"{network_type}_network.png"), dpi=300, bbox_inches='tight')
        plt.close()
    
    def create_hierarchical_network(self, network_list=None):
        """
        계층적 네트워크 모델 생성
        
        Parameters:
        -----------
        network_list : list, optional
            통합할 네트워크 유형 목록 (None이면 모든 네트워크 사용)
        """
        if network_list is None:
            network_list = list(self.networks.keys())
        
        # 계층적 네트워크 생성
        hierarchical_network = nx.DiGraph()
        
        # 각 네트워크의 노드 및 엣지 추가
        for net_type in network_list:
            if net_type not in self.networks:
                print(f"경고: '{net_type}' 네트워크를 찾을 수 없습니다.")
                continue
            
            G = self.networks[net_type]
            
            # 노드 추가
            for node, attr in G.nodes(data=True):
                # 노드 ID 충돌 방지 (네트워크 유형 접두사 추가)
                node_id = f"{net_type}:{node}"
                
                # 속성에 원본 ID와 네트워크 유형 추가
                attr_copy = attr.copy()
                attr_copy['original_id'] = node
                attr_copy['network'] = net_type
                
                hierarchical_network.add_node(node_id, **attr_copy)
            
            # 엣지 추가
            for u, v, attr in G.edges(data=True):
                u_id = f"{net_type}:{u}"
                v_id = f"{net_type}:{v}"
                
                # 속성에 네트워크 유형 추가
                attr_copy = attr.copy()
                attr_copy['network'] = net_type
                
                hierarchical_network.add_edge(u_id, v_id, **attr_copy)
        
        # 네트워크 간 연결 추가 (예: 유전자-단백질, 효소-반응 등)
        self._add_inter_network_connections(hierarchical_network)
        
        # 결과 저장
        self.networks['hierarchical'] = hierarchical_network
        return hierarchical_network
    
    def _add_inter_network_connections(self, hierarchical_network):
        """
        계층적 네트워크에서 서로 다른 네트워크 간 연결 추가
        
        Parameters:
        -----------
        hierarchical_network : NetworkX Graph
            계층적 네트워크 객체
        """
        # 1. 대사 네트워크의 효소와 PPI 네트워크의 단백질 연결
        if 'metabolic' in self.networks and 'ppi' in self.networks:
            metabolic_enzymes = [f"metabolic:{node}" for node, attr in self.networks['metabolic'].nodes(data=True) 
                                if attr.get('type') == 'enzyme']
            
            ppi_proteins = [f"ppi:{node}" for node in self.networks['ppi'].nodes()]
            
            # 효소 이름과 단백질 이름 매칭 (간단한 예시)
            for enzyme in metabolic_enzymes:
                enzyme_name = enzyme.split(':')[1]
                
                for protein in ppi_proteins:
                    protein_name = protein.split(':')[1]
                    
                    # 이름이 일치하면 연결 (실제 구현시 더 정교한 매핑 필요)
                    if enzyme_name == protein_name:
                        hierarchical_network.add_edge(enzyme, protein, type='identity', network='inter')
        
        # 2. 유전자 조절 네트워크의 유전자와 대사 네트워크의 효소 연결
        if 'gene_regulatory' in self.networks and 'metabolic' in self.networks:
            grn_genes = [f"gene_regulatory:{node}" for node, attr in self.networks['gene_regulatory'].nodes(data=True) 
                        if attr.get('type') == 'gene']
            
            metabolic_enzymes = [f"metabolic:{node}" for node, attr in self.networks['metabolic'].nodes(data=True) 
                                if attr.get('type') == 'enzyme']
            
            # 유전자-효소 매핑 (간단한 예시)
            for gene in grn_genes:
                gene_name = gene.split(':')[1]
                
                for enzyme in metabolic_enzymes:
                    enzyme_name = enzyme.split(':')[1]
                    
                    # 이름이 일치하면 연결 (실제 구현시 더 정교한 매핑 필요)
                    if gene_name == enzyme_name:
                        hierarchical_network.add_edge(gene, enzyme, type='encodes', network='inter')
        
        # 3. 신호 전달 네트워크의 전사인자와 유전자 조절 네트워크의 전사인자 연결
        if 'signaling' in self.networks and 'gene_regulatory' in self.networks:
            signaling_tfs = [f"signaling:{node}" for node, attr in self.networks['signaling'].nodes(data=True) 
                            if attr.get('type') == 'transcription_factor']
            
            grn_tfs = [f"gene_regulatory:{node}" for node, attr in self.networks['gene_regulatory'].nodes(data=True) 
                      if attr.get('type') == 'transcription_factor']
            
            # 전사인자 매핑
            for sig_tf in signaling_tfs:
                sig_tf_name = sig_tf.split(':')[1]
                
                for grn_tf in grn_tfs:
                    grn_tf_name = grn_tf.split(':')[1]
                    
                    # 이름이 일치하면 연결
                    if sig_tf_name == grn_tf_name:
                        hierarchical_network.add_edge(sig_tf, grn_tf, type='identity', network='inter')
        
        if 'drug_target' in self.networks and 'gene_regulatory' in self.networks:
            drug_genes = [f"drug_target:{n}" for n, d in self.networks['drug_target'].nodes(data=True) if d.get('type') == 'gene']
            grn_genes = [f"gene_regulatory:{n}" for n in self.networks['gene_regulatory'].nodes()]
            
            for dg in drug_genes:
                gene_id = dg.split(':')[1]
                for grn in grn_genes:
                    if grn.split(':')[1] == gene_id:
                        hierarchical_network.add_edge(dg, grn, type='target_overlap', network='inter')
    
    def visualize_hierarchical_network(self, layout='spring', highlight_inter_edges=True):
        """
        계층적 네트워크 시각화
        
        Parameters:
        -----------
        layout : str, optional
            레이아웃 알고리즘 ('spring', 'circular', 'kamada_kawai', 'spectral')
        highlight_inter_edges : bool, optional
            네트워크 간 연결 강조 여부
        """
        if 'hierarchical' not in self.networks:
            raise ValueError("먼저 계층적 네트워크를 생성하세요.")
        
        G = self.networks['hierarchical']
        
        # 네트워크 유형별 색상 매핑
        network_types = sorted({data['network'] for _, data in G.nodes(data=True)})
        network_colors = {nt: plt.cm.tab10(i) for i, nt in enumerate(network_types)}
        
        # 레이아웃 계산
        if layout == 'spring':
            pos = nx.spring_layout(G, seed=42)
        elif layout == 'circular':
            pos = nx.circular_layout(G)
        elif layout == 'kamada_kawai':
            pos = nx.kamada_kawai_layout(G)
        elif layout == 'spectral':
            pos = nx.spectral_layout(G)
        else:
            pos = nx.spring_layout(G, seed=42)  # 기본값
        
        plt.figure(figsize=(16, 14))
        
        # 네트워크별로 노드 그리기
        for net_type in network_types:
            nodes = [n for n, d in G.nodes(data=True) if d.get('network') == net_type]
            nx.draw_networkx_nodes(
                G, pos, 
                nodelist=nodes, 
                node_color=network_colors[net_type], 
                node_size=300, 
                alpha=0.7,
                label=net_type
            )
        
        # 네트워크 내부 엣지 그리기
        intra_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('network') != 'inter']
        nx.draw_networkx_edges(
            G, pos, 
            edgelist=intra_edges, 
            width=1.0, 
            alpha=0.5, 
            edge_color='gray'
        )
        
        # 네트워크 간 엣지 강조
        if highlight_inter_edges:
            inter_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('network') == 'inter']
            nx.draw_networkx_edges(
                G, pos, 
                edgelist=inter_edges, 
                width=2.0, 
                alpha=0.8, 
                edge_color='red',
                style='dashed'
            )
        
        # 범례 추가
        plt.legend()
        
        plt.title("Hierarchical Network Model")
        plt.axis('off')
        plt.tight_layout()
        
        # 결과 저장
        plt.savefig(os.path.join(self.results_dir, "hierarchical_network.png"), dpi=300, bbox_inches='tight')
        plt.close()
    
    def apply_perturbation(self, node_id, perturbation_type='knockout', network_type=None):
        """
        네트워크에 섭동 적용 (노드 활성화/억제, 엣지 삭제 등)
        
        Parameters:
        -----------
        node_id : str
            섭동을 적용할 노드 ID
        perturbation_type : str, optional
            섭동 유형 ('knockout', 'activation', 'inhibition')
        network_type : str, optional
            대상 네트워크 유형 (None이면 계층적 네트워크 사용)
        """
        if network_type is None and 'hierarchical' in self.networks:
            network_type = 'hierarchical'
        
        if network_type not in self.networks:
            raise ValueError(f"'{network_type}' 네트워크를 찾을 수 없습니다.")
        
        G = self.networks[network_type].copy()
        
        # 계층적 네트워크인 경우 네트워크 유형 접두사 처리
        if network_type == 'hierarchical' and ':' not in node_id:
            # 모든 네트워크 유형에서 노드 검색
            found_nodes = [n for n in G.nodes if n.endswith(f":{node_id}")]
            if not found_nodes:
                raise ValueError(f"노드 '{node_id}'를 찾을 수 없습니다.")
            
            perturbed_nodes = found_nodes
        else:
            # 단일 노드 처리
            if node_id not in G.nodes:
                raise ValueError(f"노드 '{node_id}'를 찾을 수 없습니다.")
            
            perturbed_nodes = [node_id]
        
        # 섭동 적용
        for node in perturbed_nodes:
            if perturbation_type == 'knockout':
                # 노드 및 연결된 엣지 제거
                G.remove_node(node)
                print(f"노드 '{node}' 제거됨")
            
            elif perturbation_type == 'activation':
                # 노드 활성화 상태 설정
                G.nodes[node]['state'] = 'active'
                G.nodes[node]['activity'] = 1.0
                print(f"노드 '{node}' 활성화됨")
            
            elif perturbation_type == 'inhibition':
                # 노드 억제 상태 설정
                G.nodes[node]['state'] = 'inhibited'
                G.nodes[node]['activity'] = 0.0
                print(f"노드 '{node}' 억제됨")
            
            else:
                raise ValueError(f"지원되지 않는 섭동 유형: {perturbation_type}")
        
        # 섭동 결과 저장
        perturbed_network_name = f"{network_type}_perturbed_{perturbation_type}_{node_id.replace(':', '_')}"
        self.networks[perturbed_network_name] = G
        
        return G
    
    def compare_networks(self, network1_type, network2_type, metric='node_overlap'):
        """
        두 네트워크 비교
        
        Parameters:
        -----------
        network1_type : str
            첫 번째 네트워크 유형
        network2_type : str
            두 번째 네트워크 유형
        metric : str, optional
            비교 지표 ('node_overlap', 'edge_overlap', 'jaccard', 'graph_edit_distance')
        """
        if network1_type not in self.networks:
            raise ValueError(f"'{network1_type}' 네트워크를 찾을 수 없습니다.")
        
        if network2_type not in self.networks:
            raise ValueError(f"'{network2_type}' 네트워크를 찾을 수 없습니다.")
        
        G1 = self.networks[network1_type]
        G2 = self.networks[network2_type]
        
        # 비교 지표 계산
        results = {}
        
        if metric == 'node_overlap' or metric == 'all':
            # 노드 집합 비교
            nodes1 = set(G1.nodes())
            nodes2 = set(G2.nodes())
            
            common_nodes = nodes1.intersection(nodes2)
            only_in_1 = nodes1 - nodes2
            only_in_2 = nodes2 - nodes1
            
            results['node_overlap'] = {
                'common_nodes': len(common_nodes),
                'only_in_network1': len(only_in_1),
                'only_in_network2': len(only_in_2),
                'jaccard_nodes': len(common_nodes) / len(nodes1.union(nodes2)) if nodes1 or nodes2 else 0
            }
        
        if metric == 'edge_overlap' or metric == 'all':
            # 엣지 집합 비교
            edges1 = set(G1.edges())
            edges2 = set(G2.edges())
            
            common_edges = edges1.intersection(edges2)
            only_in_1 = edges1 - edges2
            only_in_2 = edges2 - edges1
            
            results['edge_overlap'] = {
                'common_edges': len(common_edges),
                'only_in_network1': len(only_in_1),
                'only_in_network2': len(only_in_2),
                'jaccard_edges': len(common_edges) / len(edges1.union(edges2)) if edges1 or edges2 else 0
            }
        
        if metric == 'jaccard' or metric == 'all':
            # 노드 자카드 유사도
            nodes1 = set(G1.nodes())
            nodes2 = set(G2.nodes())
            jaccard_nodes = len(nodes1.intersection(nodes2)) / len(nodes1.union(nodes2)) if nodes1 or nodes2 else 0
            
            # 엣지 자카드 유사도
            edges1 = set(G1.edges())
            edges2 = set(G2.edges())
            jaccard_edges = len(edges1.intersection(edges2)) / len(edges1.union(edges2)) if edges1 or edges2 else 0
            
            results['jaccard'] = {
                'jaccard_nodes': jaccard_nodes,
                'jaccard_edges': jaccard_edges
            }
        
        if metric == 'graph_edit_distance' or metric == 'all':
            # 그래프 편집 거리 계산 (계산에 시간이 오래 걸릴 수 있음)
            try:
                # 약간의 근사치 알고리즘 사용 (더 빠르지만 정확하지 않을 수 있음)
                ged = nx.graph_edit_distance(G1, G2, timeout=10)  # 10초 타임아웃
                results['graph_edit_distance'] = ged
            except:
                print("그래프 편집 거리 계산 시간 초과 또는 오류 발생")
                results['graph_edit_distance'] = "계산 불가"
        
        return results
    
    #----------------------------------------------------
    # 3. 네트워크 모듈 분석
    #----------------------------------------------------
    
    def identify_modules(self, network_type, algorithm='louvain', resolution=1.0, k=None):
        """
        네트워크 모듈 식별
        
        Parameters:
        -----------
        network_type : str
            모듈을 식별할 네트워크 유형
        algorithm : str, optional
            커뮤니티 검출 알고리즘 ('louvain', 'kmeans', 'spectral')
        resolution : float, optional
            모듈 해상도 파라미터 (Louvain 알고리즘용)
        k : int, optional
            모듈 수 (k-means, spectral 알고리즘용)
        """
        if network_type not in self.networks:
            raise ValueError(f"'{network_type}' 네트워크를 찾을 수 없습니다.")
        
        G = self.networks[network_type]
        
        # 알고리즘에 따른 모듈 식별
        if algorithm == 'louvain':
            # Louvain 커뮤니티 검출 알고리즘
            if nx.is_directed(G):
                # 방향성 그래프는 무방향 그래프로 변환
                undirected_G = G.to_undirected()
                partition = community_louvain.best_partition(undirected_G, resolution=resolution)
            else:
                partition = community_louvain.best_partition(G, resolution=resolution)
            
            # 결과 포맷팅
            modules = {}
            for node, module_id in partition.items():
                if module_id not in modules:
                    modules[module_id] = []
                modules[module_id].append(node)
            
            # 모듈 정보를 네트워크에 속성으로 추가
            nx.set_node_attributes(G, partition, 'module')
        
        elif algorithm == 'kmeans':
            if k is None:
                # 기본 k 값 설정: 노드 수의 제곱근
                k = int(np.sqrt(G.number_of_nodes()))
            
            # 네트워크를 행렬로 변환
            adj_matrix = nx.to_numpy_array(G)
            
            # K-means 클러스터링
            kmeans = KMeans(n_clusters=k, random_state=42)
            node_list = list(G.nodes())
            
            if adj_matrix.shape[0] > 1:  # 행렬이 비어있지 않은 경우
                cluster_labels = kmeans.fit_predict(adj_matrix)
                
                # 결과 포맷팅
                modules = {}
                for i, label in enumerate(cluster_labels):
                    if label not in modules:
                        modules[label] = []
                    modules[label].append(node_list[i])
                
                # 모듈 정보를 네트워크에 속성으로 추가
                module_dict = {node_list[i]: label for i, label in enumerate(cluster_labels)}
                nx.set_node_attributes(G, module_dict, 'module')
            else:
                modules = {}
        
        elif algorithm == 'spectral':
            if k is None:
                # 기본 k 값 설정
                k = min(10, G.number_of_nodes())
            
            # 네트워크를 행렬로 변환
            adj_matrix = nx.to_numpy_array(G)
            
            # 행렬이 비어있지 않고 충분한 크기인 경우
            if adj_matrix.shape[0] > k:
                # Spectral 클러스터링
                spectral = SpectralClustering(n_clusters=k, random_state=42, affinity='precomputed')
                node_list = list(G.nodes())
                
                # 인접 행렬이 희소행렬인 경우 밀집행렬로 변환
                if isinstance(adj_matrix, np.ndarray) and adj_matrix.size > 0:
                    cluster_labels = spectral.fit_predict(adj_matrix)
                    
                    # 결과 포맷팅
                    modules = {}
                    for i, label in enumerate(cluster_labels):
                        if label not in modules:
                            modules[label] = []
                        modules[label].append(node_list[i])
                    
                    # 모듈 정보를 네트워크에 속성으로 추가
                    module_dict = {node_list[i]: label for i, label in enumerate(cluster_labels)}
                    nx.set_node_attributes(G, module_dict, 'module')
                else:
                    modules = {}
            else:
                print(f"경고: 노드 수({G.number_of_nodes()})가 k({k})보다 작거나 같아 Spectral 클러스터링을 적용할 수 없습니다.")
                modules = {}
        
        else:
            raise ValueError(f"지원되지 않는 알고리즘: {algorithm}")
        
        # 모듈 결과 저장
        self.modules[network_type] = {
            'algorithm': algorithm,
            'modules': modules,
            'node_to_module': {node: module_id for module_id, nodes in modules.items() for node in nodes} if modules else {}
        }
        
        return modules
    
    def visualize_modules(self, network_type, layout='spring'):
        """
        네트워크 모듈 시각화
        
        Parameters:
        -----------
        network_type : str
            시각화할 네트워크 유형
        layout : str, optional
            레이아웃 알고리즘 ('spring', 'circular', 'kamada_kawai', 'spectral')
        """
        if network_type not in self.networks:
            raise ValueError(f"'{network_type}' 네트워크를 찾을 수 없습니다.")
        
        if network_type not in self.modules:
            raise ValueError(f"'{network_type}' 네트워크의 모듈 분석 결과를 찾을 수 없습니다.")
        
        G = self.networks[network_type]
        modules_info = self.modules[network_type]
        
        # 모듈 정보 확인
        if not modules_info['modules']:
            raise ValueError("모듈 정보가 비어 있습니다.")
        
        # 레이아웃 계산
        if layout == 'spring':
            pos = nx.spring_layout(G, seed=42)
        elif layout == 'circular':
            pos = nx.circular_layout(G)
        elif layout == 'kamada_kawai':
            pos = nx.kamada_kawai_layout(G)
        elif layout == 'spectral':
            pos = nx.spectral_layout(G)
        else:
            pos = nx.spring_layout(G, seed=42)  # 기본값
        
        plt.figure(figsize=(14, 12))
        
        # 모듈별 노드 색상 설정
        node_colors = []
        module_ids = sorted(modules_info['modules'].keys())
        color_map = {module_id: plt.cm.tab20(i % 20) for i, module_id in enumerate(module_ids)}
        
        for node in G.nodes():
            module_id = G.nodes[node].get('module', -1)
            node_colors.append(color_map.get(module_id, 'gray'))
        
        # 노드 그리기
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=300, alpha=0.8)
        
        # 엣지 그리기
        nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.5, edge_color='gray')
        
        # 노드 레이블 그리기
        nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold')
        
        # 범례 추가
        legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map[module_id], 
                                     label=f'Module {module_id}', markersize=10) 
                         for module_id in module_ids]
        plt.legend(handles=legend_elements, loc='upper right')
        
        plt.title(f"Modules in {network_type.capitalize()} Network")
        plt.axis('off')
        plt.tight_layout()
        
        # 결과 저장
        plt.savefig(os.path.join(self.results_dir, f"{network_type}_modules.png"), dpi=300, bbox_inches='tight')
        plt.close()
    
    def analyze_module_functions(self, network_type, module_id=None, gene_annotation=None):
        """
        모듈 기능 분석
        
        Parameters:
        -----------
        network_type : str
            분석할 네트워크 유형
        module_id : int or str, optional
            분석할 특정 모듈 ID (None이면 모든 모듈 분석)
        gene_annotation : dict or pd.DataFrame, optional
            유전자/단백질 주석 정보
        """
        if network_type not in self.modules:
            raise ValueError(f"'{network_type}' 네트워크의 모듈 분석 결과를 찾을 수 없습니다.")
        
        modules_info = self.modules[network_type]
        modules = modules_info['modules']
        
        # 분석할 모듈 선택
        if module_id is not None:
            if module_id not in modules:
                raise ValueError(f"모듈 ID '{module_id}'를 찾을 수 없습니다.")
            
            modules_to_analyze = {module_id: modules[module_id]}
        else:
            modules_to_analyze = modules
        
        # 유전자 주석 확인
        if gene_annotation is None:
            # 예시 주석 데이터 생성 (실제 구현시 제거 필요)
            gene_annotation = {
                'Gene1': {'function': 'Transcription factor', 'pathway': 'Cell cycle', 'GO': ['GO:0006351']},
                'Gene2': {'function': 'Kinase', 'pathway': 'MAPK signaling', 'GO': ['GO:0004672']},
                'Gene3': {'function': 'Transporter', 'pathway': 'ABC transporters', 'GO': ['GO:0005215']},
                'Protein1': {'function': 'Enzyme', 'pathway': 'Glycolysis', 'GO': ['GO:0016491']},
                'Protein2': {'function': 'Receptor', 'pathway': 'GPCR signaling', 'GO': ['GO:0004930']},
                'Protein3': {'function': 'Chaperone', 'pathway': 'Protein folding', 'GO': ['GO:0006457']},
                'HK': {'function': 'Hexokinase', 'pathway': 'Glycolysis', 'GO': ['GO:0004396']},
                'PGI': {'function': 'Phosphoglucose isomerase', 'pathway': 'Glycolysis', 'GO': ['GO:0004347']},
                'TF1': {'function': 'Transcription factor', 'pathway': 'Gene regulation', 'GO': ['GO:0003700']},
                'TF2': {'function': 'Transcription factor', 'pathway': 'Gene regulation', 'GO': ['GO:0003700']},
                'Kinase1': {'function': 'Kinase', 'pathway': 'Signal transduction', 'GO': ['GO:0004672']},
                'Kinase2': {'function': 'Kinase', 'pathway': 'Signal transduction', 'GO': ['GO:0004672']},
            }
        
        # 모듈 기능 분석 결과
        module_functions = {}
        
        for mid, nodes in modules_to_analyze.items():
            # 모듈 내 노드 분석
            functions = {}
            pathways = {}
            go_terms = {}
            
            for node in nodes:
                # 네트워크 유형 접두사 제거 (계층적 네트워크의 경우)
                if ':' in node:
                    node_name = node.split(':')[1]
                else:
                    node_name = node
                
                # 유전자/단백질 주석 검색
                annotation = gene_annotation.get(node_name, {})
                
                # 기능 집계
                function = annotation.get('function', 'Unknown')
                if function in functions:
                    functions[function] += 1
                else:
                    functions[function] = 1
                
                # 경로 집계
                pathway = annotation.get('pathway', 'Unknown')
                if pathway in pathways:
                    pathways[pathway] += 1
                else:
                    pathways[pathway] = 1
                
                # GO 용어 집계
                for go in annotation.get('GO', []):
                    if go in go_terms:
                        go_terms[go] += 1
                    else:
                        go_terms[go] = 1
            
            # 가장 빈번한 기능/경로/GO 용어 추출
            top_functions = sorted(functions.items(), key=lambda x: x[1], reverse=True)
            top_pathways = sorted(pathways.items(), key=lambda x: x[1], reverse=True)
            top_go_terms = sorted(go_terms.items(), key=lambda x: x[1], reverse=True)
            
            module_functions[mid] = {
                'node_count': len(nodes),
                'top_functions': top_functions[:5],  # 상위 5개만
                'top_pathways': top_pathways[:5],    # 상위 5개만
                'top_go_terms': top_go_terms[:5],    # 상위 5개만
            }
        
        # 결과 저장
        if module_id is not None:
            self.modules[network_type]['module_functions'] = {module_id: module_functions[module_id]}
        else:
            self.modules[network_type]['module_functions'] = module_functions
        
        return module_functions
    
    def export_module_analysis(self, network_type, output_file=None):
        """
        모듈 분석 결과 내보내기
        
        Parameters:
        -----------
        network_type : str
            내보낼 네트워크 유형
        output_file : str, optional
            출력 파일 경로 (None이면 기본 경로 사용)
        """
        if network_type not in self.modules:
            raise ValueError(f"'{network_type}' 네트워크의 모듈 분석 결과를 찾을 수 없습니다.")
        
        modules_info = self.modules[network_type]
        
        # 기본 출력 파일 경로
        if output_file is None:
            output_file = os.path.join(self.results_dir, f"{network_type}_module_analysis.json")
        
        # 결과를 JSON 형식으로 내보내기
        with open(output_file, 'w') as f:
            json.dump(modules_info, f, indent=2)
        
        print(f"모듈 분석 결과가 '{output_file}'에 저장되었습니다.")
        
        # 모듈 요약 보고서 생성
        report_file = os.path.join(self.results_dir, f"{network_type}_module_report.txt")
        
        with open(report_file, 'w') as f:
            f.write(f"# {network_type.capitalize()} 네트워크 모듈 분석 보고서\n\n")
            f.write(f"알고리즘: {modules_info['algorithm']}\n")
            f.write(f"모듈 수: {len(modules_info['modules'])}\n\n")
            
            # 각 모듈에 대한 정보
            for module_id, nodes in modules_info['modules'].items():
                f.write(f"## 모듈 {module_id}\n")
                f.write(f"노드 수: {len(nodes)}\n")
                
                # 모듈 내 상위 노드 (최대 10개)
                f.write("주요 노드:\n")
                for node in nodes[:10]:
                    f.write(f"- {node}\n")
                
                if len(nodes) > 10:
                    f.write(f"... 외 {len(nodes) - 10}개 노드\n")
                
                # 모듈 기능 정보 (있는 경우)
                if 'module_functions' in modules_info and module_id in modules_info['module_functions']:
                    func_info = modules_info['module_functions'][module_id]
                    
                    f.write("\n주요 기능:\n")
                    for func, count in func_info.get('top_functions', []):
                        f.write(f"- {func}: {count}개 노드\n")
                    
                    f.write("\n주요 경로:\n")
                    for pathway, count in func_info.get('top_pathways', []):
                        f.write(f"- {pathway}: {count}개 노드\n")
                
                f.write("\n")
        
        print(f"모듈 분석 보고서가 '{report_file}'에 저장되었습니다.")
    
    #----------------------------------------------------
    # 4. 시간적 네트워크 역학
    #----------------------------------------------------
    
    def load_time_series_data(self, file_path, time_col=None, id_col=None):
        """
        시계열 데이터 로드
        
        Parameters:
        -----------
        file_path : str
            시계열 데이터 파일 경로 (CSV, TSV 등)
        time_col : str, optional
            시간 포인트 열 이름
        id_col : str, optional
            특성 ID 열 이름
        """
        # 데이터 로드
        if file_path.endswith('.csv'):
            df = pd.read_csv(file_path)
        elif file_path.endswith('.tsv') or file_path.endswith('.txt'):
            df = pd.read_csv(file_path, sep='\t')
        elif file_path.endswith('.xlsx'):
            df = pd.read_excel(file_path)
        else:
            raise ValueError("지원되지 않는 파일 형식입니다.")
        
        # 예시 데이터 생성 (실제 구현시 제거 필요)
        if 'example' in file_path.lower():
            # 예시 시계열 데이터 생성
            np.random.seed(42)
            times = np.arange(0, 24, 2)  # 0, 2, 4, ..., 22 시간
            genes = [f"Gene{i}" for i in range(1, 21)]  # Gene1, Gene2, ..., Gene20
            
            data = []
            for gene in genes:
                # 기본 발현 패턴 생성 (사인파 또는 임의 패턴)
                if np.random.rand() < 0.5:
                    # 사인파 패턴
                    period = np.random.uniform(6, 24)  # 주기: 6~24시간
                    phase = np.random.uniform(0, 2*np.pi)  # 위상
                    amplitude = np.random.uniform(0.5, 2.0)  # 진폭
                    base = np.random.uniform(0.5, 1.5)  # 기준값
                    
                    expression = base + amplitude * np.sin(2*np.pi*times/period + phase)
                else:
                    # 임의 패턴
                    expression = np.random.normal(1.0, 0.3, size=len(times))
                    # 연속적인 추세 추가
                    for i in range(1, len(expression)):
                        expression[i] = 0.7 * expression[i] + 0.3 * expression[i-1]
                
                # 약간의 노이즈 추가
                expression += np.random.normal(0, 0.1, size=len(times))
                
                # 음수 방지
                expression = np.maximum(0, expression)
                
                # 데이터 추가
                for i, t in enumerate(times):
                    data.append({'Time': t, 'GeneID': gene, 'Expression': expression[i]})
            
            df = pd.DataFrame(data)
            time_col = 'Time'
            id_col = 'GeneID'
        
        # 시간 열 및 ID 열 처리
        if time_col is not None and time_col in df.columns:
            time_points = sorted(df[time_col].unique())
        else:
            # 시간 열이 없으면 열 이름을 시간 포인트로 가정
            numeric_cols = [col for col in df.columns if col != id_col and pd.api.types.is_numeric_dtype(df[col])]
            time_points = sorted([float(col) for col in numeric_cols if col.replace('.', '').isdigit()])
        
        # 결과 포맷팅
        if time_col is not None and id_col is not None:
            # 세로 형식 데이터 (long format)
            # 시간을 기준으로 피벗
            pivoted_df = df.pivot(index=id_col, columns=time_col, values='Expression')
            time_series_data = pivoted_df
        elif id_col is not None:
            # 가로 형식 데이터 (wide format)
            time_series_data = df.set_index(id_col)
        else:
            # 모든 열이 시간 포인트로 가정
            time_series_data = df
        
        # 시계열 데이터 저장
        self.dynamics['time_series_data'] = time_series_data
        self.dynamics['time_points'] = time_points
        
        return time_series_data
    
    def create_dynamic_network(self, network_type, time_series_data=None, threshold=0.5):
        """
        시간적 네트워크 생성
        
        Parameters:
        -----------
        network_type : str
            기반 네트워크 유형
        time_series_data : pd.DataFrame, optional
            시계열 발현 데이터 (None이면 저장된 데이터 사용)
        threshold : float, optional
            활성화 임계값
        """
        if network_type not in self.networks:
            raise ValueError(f"'{network_type}' 네트워크를 찾을 수 없습니다.")
        
        G = self.networks[network_type]
        
        # 시계열 데이터 확인
        if time_series_data is None:
            if 'time_series_data' not in self.dynamics:
                raise ValueError("시계열 데이터를 먼저 로드하세요.")
            
            time_series_data = self.dynamics['time_series_data']
            time_points = self.dynamics['time_points']
        else:
            # 시간 포인트 추출
            time_points = time_series_data.columns
        
        # 동적 네트워크 생성
        dynamic_networks = {}
        
        for t in time_points:
            # 현재 시간 포인트의 발현값
            expression_t = time_series_data[t]
            
            # 정규화 (필요시)
            if expression_t.max() > 1.0:
                expression_t = expression_t / expression_t.max()
            
            # t 시점의 네트워크 복사
            G_t = G.copy()
            
            # 노드 활성화 상태 업데이트
            for node in G_t.nodes():
                # 네트워크 유형 접두사 제거 (계층적 네트워크의 경우)
                if ':' in node:
                    node_name = node.split(':')[1]
                else:
                    node_name = node
                
                # 발현 데이터에서 노드 검색
                if node_name in expression_t.index:
                    # 활성화 상태 설정
                    G_t.nodes[node]['activity'] = float(expression_t[node_name])
                    G_t.nodes[node]['state'] = 'active' if expression_t[node_name] >= threshold else 'inactive'
                else:
                    # 데이터가 없는 경우 기본값
                    G_t.nodes[node]['activity'] = 0.0
                    G_t.nodes[node]['state'] = 'unknown'
            
            # 엣지 가중치 업데이트 (활성화된 노드 간의 연결 강화)
            for u, v, data in G_t.edges(data=True):
                u_activity = G_t.nodes[u].get('activity', 0.0)
                v_activity = G_t.nodes[v].get('activity', 0.0)
                
                # 가중치 업데이트 (노드 활성도의 곱)
                G_t[u][v]['weight'] = u_activity * v_activity
            
            dynamic_networks[t] = G_t
        
        # 동적 네트워크 저장
        self.dynamics['dynamic_networks'] = dynamic_networks
        self.dynamics['network_type'] = network_type
        
        return dynamic_networks
    
    def simulate_dynamic_network(self, initial_state=None, steps=10, model='discrete'):
        """
        동적 네트워크 시뮬레이션
        
        Parameters:
        -----------
        initial_state : dict, optional
            초기 노드 상태 (None이면 무작위 초기화)
        steps : int, optional
            시뮬레이션 단계 수
        model : str, optional
            시뮬레이션 모델 ('discrete', 'continuous')
        """
        if 'network_type' not in self.dynamics:
            raise ValueError("먼저 동적 네트워크를 생성하세요.")
        
        network_type = self.dynamics['network_type']
        G = self.networks[network_type]
        
        # 초기 상태 설정
        if initial_state is None:
            # 무작위 초기화
            initial_state = {}
            for node in G.nodes():
                initial_state[node] = np.random.rand()
        
        # 노드 상태 업데이트
        node_states = {0: initial_state}
        
        if model == 'discrete':
            # 이산 시간 모델
            current_state = initial_state.copy()
            
            for t in range(1, steps + 1):
                next_state = {}
                
                for node in G.nodes():
                    # 각 노드의 상태 계산
                    incoming_nodes = list(G.predecessors(node))
                    
                    if not incoming_nodes:
                        # 입력이 없는 노드는 상태 유지
                        next_state[node] = current_state[node]
                    else:
                        # 입력 노드 상태 가중 합
                        weighted_sum = 0.0
                        total_weight = 0.0
                        
                        for u in incoming_nodes:
                            weight = G[u][node].get('weight', 1.0)
                            weighted_sum += current_state[u] * weight
                            total_weight += weight
                        
                        if total_weight > 0:
                            avg_input = weighted_sum / total_weight
                        else:
                            avg_input = 0.0
                        
                        # 시그모이드 활성화 함수
                        next_state[node] = 1.0 / (1.0 + np.exp(-4 * (avg_input - 0.5)))
                
                # 다음 상태 저장
                node_states[t] = next_state
                current_state = next_state.copy()
        
        elif model == 'continuous':
            # 연속 시간 모델 (간단한 ODE 시스템)
            
            # 미분 방정식 정의
            def network_dynamics(t, y, G, node_list):
                dy = np.zeros_like(y)
                
                for i, node in enumerate(node_list):
                    # 각 노드의 변화율 계산
                    incoming_nodes = list(G.predecessors(node))
                    
                    if not incoming_nodes:
                        # 입력이 없는 노드는 감쇠
                        dy[i] = -0.1 * y[i]
                    else:
                        # 입력 노드의 영향
                        input_effect = 0.0
                        
                        for u in incoming_nodes:
                            u_idx = node_list.index(u)
                            weight = G[u][node].get('weight', 1.0)
                            
                            # 엣지 유형에 따른
                            if 'effect' in G[u][node]:
                                if G[u][node]['effect'] == 'activation':
                                    input_effect += weight * y[u_idx]
                                elif G[u][node]['effect'] == 'inhibition' or G[u][node]['effect'] == 'repression':
                                    input_effect -= weight * y[u_idx]
                                else:
                                    input_effect += weight * y[u_idx]
                            else:
                                input_effect += weight * y[u_idx]
                        
                        # 변화율 계산 (성장 및 감쇠 항 포함)
                        dy[i] = 0.1 * input_effect - 0.1 * y[i]
                
                return dy
            
            # 노드 리스트 및 초기 상태 벡터
            node_list = list(G.nodes())
            y0 = np.array([initial_state[node] for node in node_list])
            
            # 시간 포인트
            t_span = (0, steps)
            t_eval = np.linspace(0, steps, steps + 1)
            
            # ODE 해석
            sol = solve_ivp(
                lambda t, y: network_dynamics(t, y, G, node_list),
                t_span, y0, t_eval=t_eval, method='RK45'
            )
            
            # 결과 변환
            for i, t in enumerate(sol.t):
                node_states[int(t)] = {node: max(0, min(1, sol.y[j][i])) for j, node in enumerate(node_list)}
        
        else:
            raise ValueError(f"지원되지 않는 모델: {model}")
        
        # 시뮬레이션 결과 저장
        self.dynamics['simulation'] = {
            'model': model,
            'steps': steps,
            'node_states': node_states
        }
        
        return node_states
    
    def visualize_dynamic_network(self, time_points=None, layout='spring'):
        """
        동적 네트워크 시각화
        
        Parameters:
        -----------
        time_points : list, optional
            시각화할 시간 포인트 (None이면 모든 시간 포인트)
        layout : str, optional
            레이아웃 알고리즘 ('spring', 'circular', 'kamada_kawai', 'spectral')
        """
        if 'dynamic_networks' not in self.dynamics:
            raise ValueError("먼저 동적 네트워크를 생성하세요.")
        
        dynamic_networks = self.dynamics['dynamic_networks']
        
        # 시각화할 시간 포인트 선택
        if time_points is None:
            time_points = sorted(dynamic_networks.keys())
        else:
            # 존재하는 시간 포인트만 선택
            time_points = sorted([t for t in time_points if t in dynamic_networks])
        
        if not time_points:
            raise ValueError("시각화할 시간 포인트가 없습니다.")
        
        # 모든 시간 포인트에 대해 동일한 레이아웃 사용
        G_first = dynamic_networks[time_points[0]]
        
        if layout == 'spring':
            pos = nx.spring_layout(G_first, seed=42)
        elif layout == 'circular':
            pos = nx.circular_layout(G_first)
        elif layout == 'kamada_kawai':
            pos = nx.kamada_kawai_layout(G_first)
        elif layout == 'spectral':
            pos = nx.spectral_layout(G_first)
        else:
            pos = nx.spring_layout(G_first, seed=42)  # 기본값
        
        # 각 시간 포인트별 네트워크 시각화
        for t in time_points:
            G_t = dynamic_networks[t]
            
            plt.figure(figsize=(12, 10))
            
            # 노드 활성도에 따른 색상 및 크기 설정
            node_colors = []
            node_sizes = []
            
            for node in G_t.nodes():
                activity = G_t.nodes[node].get('activity', 0.0)
                
                # 활성도에 따른 색상 (0: 파란색, 1: 빨간색)
                node_colors.append(plt.cm.coolwarm(activity))
                
                # 활성도에 따른 크기
                node_sizes.append(300 + 200 * activity)
            
            # 노드 그리기
            nx.draw_networkx_nodes(G_t, pos, node_color=node_colors, node_size=node_sizes, alpha=0.8)
            
            # 엣지 그리기 (가중치에 따른 두께)
            edges = [(u, v) for u, v, d in G_t.edges(data=True)]
            edge_weights = [d.get('weight', 1.0) for _, _, d in G_t.edges(data=True)]
            
            nx.draw_networkx_edges(
                G_t, pos, 
                edgelist=edges, 
                width=[max(0.5, 3 * w) for w in edge_weights], 
                alpha=0.6, 
                edge_color='gray'
            )
            
            # 노드 레이블 그리기
            nx.draw_networkx_labels(G_t, pos, font_size=8, font_weight='bold')
            
            plt.title(f"Dynamic Network at Time {t}")
            plt.axis('off')
            plt.tight_layout()
            
            # 결과 저장
            plt.savefig(os.path.join(self.results_dir, f"dynamic_network_t{t}.png"), dpi=300, bbox_inches='tight')
            plt.close()
        
        # 시간에 따른 노드 활성도 변화 플롯
        plt.figure(figsize=(14, 8))
        
        # 몇 개의 주요 노드 선택
        sample_nodes = list(G_first.nodes())[:10]  # 처음 10개 노드
        
        for node in sample_nodes:
            activities = [dynamic_networks[t].nodes[node].get('activity', 0.0) for t in time_points]
            plt.plot(time_points, activities, marker='o', label=node)
        
        plt.xlabel('Time')
        plt.ylabel('Node Activity')
        plt.title('Node Activity Changes Over Time')
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.7)
        
        # 결과 저장
        plt.savefig(os.path.join(self.results_dir, "node_activity_time_series.png"), dpi=300, bbox_inches='tight')
        plt.close()
    
    def analyze_dynamic_patterns(self, clustering=True, n_clusters=3):
        """
        시간적 패턴 분석
        
        Parameters:
        -----------
        clustering : bool, optional
            시간적 패턴 클러스터링 수행 여부
        n_clusters : int, optional
            클러스터 수
        """
        if 'dynamic_networks' not in self.dynamics:
            raise ValueError("먼저 동적 네트워크를 생성하세요.")
        
        dynamic_networks = self.dynamics['dynamic_networks']
        time_points = sorted(dynamic_networks.keys())
        
        # 활성도 시계열 데이터 생성
        nodes = list(dynamic_networks[time_points[0]].nodes())
        activity_data = {}
        
        for node in nodes:
            activity_data[node] = [dynamic_networks[t].nodes[node].get('activity', 0.0) for t in time_points]
        
        # 시계열 데이터 프레임 생성
        activity_df = pd.DataFrame(activity_data, index=time_points).T
        
        # 패턴 분석 결과
        pattern_analysis = {
            'activity_df': activity_df,
            'time_points': time_points
        }
        
        # 클러스터링 수행
        if clustering:
            # 클러스터 수가 데이터 수보다 작은지 확인
            n_clusters = min(n_clusters, len(nodes))
            
            # K-means 클러스터링
            kmeans = KMeans(n_clusters=n_clusters, random_state=42)
            clusters = kmeans.fit_predict(activity_df.values)
            
            # 클러스터 결과 추가
            pattern_analysis['clusters'] = pd.Series(clusters, index=activity_df.index)
            pattern_analysis['n_clusters'] = n_clusters
            pattern_analysis['cluster_centers'] = kmeans.cluster_centers_
            
            # 클러스터별 평균 활성도 패턴
            cluster_means = {}
            for i in range(n_clusters):
                cluster_nodes = activity_df.index[clusters == i]
                cluster_means[i] = activity_df.loc[cluster_nodes].mean()
            
            pattern_analysis['cluster_means'] = pd.DataFrame(cluster_means)
        
        # 결과 저장
        self.dynamics['pattern_analysis'] = pattern_analysis
        
        return pattern_analysis
    
    def visualize_dynamic_patterns(self):
        """시간적 패턴 시각화"""
        if 'pattern_analysis' not in self.dynamics:
            raise ValueError("먼저 시간적 패턴 분석을 수행하세요.")
        
        pattern_analysis = self.dynamics['pattern_analysis']
        activity_df = pattern_analysis['activity_df']
        time_points = pattern_analysis['time_points']
        
        # 1. 히트맵 시각화
        plt.figure(figsize=(12, 10))
        sns.heatmap(
            activity_df,
            cmap='coolwarm',
            vmin=0,
            vmax=1,
            xticklabels=time_points,
            yticklabels=False,  # 노드가 많을 경우 레이블 생략
            cbar_kws={'label': 'Activity'}
        )
        plt.title('Node Activity Patterns Over Time')
        plt.xlabel('Time')
        plt.ylabel('Nodes')
        plt.tight_layout()
        
        # 결과 저장
        plt.savefig(os.path.join(self.results_dir, "activity_heatmap.png"), dpi=300, bbox_inches='tight')
        plt.close()
        
        # 2. 클러스터링 결과 시각화 (수행한 경우)
        if 'clusters' in pattern_analysis:
            clusters = pattern_analysis['clusters']
            n_clusters = pattern_analysis['n_clusters']
            cluster_means = pattern_analysis['cluster_means']
            
            # 클러스터별 평균 활성도 패턴
            plt.figure(figsize=(12, 8))
            for i in range(n_clusters):
                plt.plot(
                    time_points, 
                    cluster_means[i], 
                    marker='o',
                    linewidth=2,
                    label=f'Cluster {i} (n={sum(clusters == i)})'
                )
            
            plt.xlabel('Time')
            plt.ylabel('Average Activity')
            plt.title('Cluster Mean Activity Patterns')
            plt.legend()
            plt.grid(True, linestyle='--', alpha=0.7)
            
            # 결과 저장
            plt.savefig(os.path.join(self.results_dir, "cluster_means.png"), dpi=300, bbox_inches='tight')
            plt.close()
            
            # 3. 개별 클러스터 히트맵
            for i in range(n_clusters):
                # 클러스터에 속한 노드 추출
                cluster_nodes = activity_df.index[clusters == i]
                
                # 클러스터 크기가 너무 크면 샘플링
                if len(cluster_nodes) > 50:
                    print(f"클러스터 {i}의 노드가 너무 많습니다. 처음 50개만 시각화합니다.")
                    cluster_nodes = cluster_nodes[:50]
                
                # 충분한 노드가 있는 경우에만 시각화
                if len(cluster_nodes) > 0:
                    plt.figure(figsize=(12, min(10, max(6, len(cluster_nodes) / 5))))
                    sns.heatmap(
                        activity_df.loc[cluster_nodes],
                        cmap='coolwarm',
                        vmin=0,
                        vmax=1,
                        xticklabels=time_points,
                        yticklabels=True,  # 클러스터 내 노드 레이블 표시
                        cbar_kws={'label': 'Activity'}
                    )
                    plt.title(f'Activity Patterns for Cluster {i} (n={sum(clusters == i)})')
                    plt.xlabel('Time')
                    plt.ylabel('Nodes')
                    plt.tight_layout()
                    
                    # 결과 저장
                    plt.savefig(os.path.join(self.results_dir, f"cluster_{i}_heatmap.png"), dpi=300, bbox_inches='tight')
                    plt.close()
            
            # 4. UMAP 또는 t-SNE를 통한 시계열 패턴 2D 시각화
            try:
                # UMAP 차원 축소
                reducer = umap.UMAP(n_components=2, random_state=42)
                embedding = reducer.fit_transform(activity_df.values)
                
                plt.figure(figsize=(10, 8))
                
                # 클러스터별 색상
                colors = plt.cm.tab10(np.array([i for i in range(n_clusters)]))
                
                for i in range(n_clusters):
                    cluster_mask = clusters == i
                    plt.scatter(
                        embedding[cluster_mask, 0],
                        embedding[cluster_mask, 1],
                        color=colors[i],
                        label=f'Cluster {i}',
                        alpha=0.7,
                        s=80
                    )
                
                plt.title('UMAP Visualization of Activity Patterns')
                plt.legend()
                plt.tight_layout()
                
                # 결과 저장
                plt.savefig(os.path.join(self.results_dir, "umap_activity_patterns.png"), dpi=300, bbox_inches='tight')
                plt.close()
            except:
                print("UMAP 시각화 중 오류 발생")
    
    #----------------------------------------------------
    # 5. 다중 스케일 모델링
    #----------------------------------------------------
    
    def create_multiscale_model(self, model_scales=['molecular', 'cellular'], integration_method='vertical'):
        """
        다중 스케일 모델 생성
        
        Parameters:
        -----------
        model_scales : list, optional
            포함할 모델 스케일 목록
        integration_method : str, optional
            통합 방법 ('vertical', 'horizontal', 'nested')
        """
        # 다중 스케일 모델 구조
        multiscale_model = {
            'scales': {},
            'integration_method': integration_method,
            'scale_connections': []
        }
        
        # 각 스케일별 모델 생성
        for scale in model_scales:
            if scale == 'molecular':
                # 분자 수준 모델 (대사 네트워크, PPI 등)
                scale_model = {
                    'networks': {},
                    'parameters': {}
                }
                
                # 네트워크 포함
                for net_type in ['metabolic', 'ppi']:
                    if net_type in self.networks:
                        scale_model['networks'][net_type] = self.networks[net_type]
                
                # 분자 수준 매개변수 설정
                scale_model['parameters'] = {
                    'diffusion_rate': 0.1,
                    'binding_constant': 1.0,
                    'reaction_rates': {'default': 0.5}
                }
                
                multiscale_model['scales']['molecular'] = scale_model
            
            elif scale == 'cellular':
                # 세포 수준 모델 (유전자 조절 네트워크, 신호 전달 등)
                scale_model = {
                    'networks': {},
                    'parameters': {}
                }
                
                # 네트워크 포함
                for net_type in ['gene_regulatory', 'signaling']:
                    if net_type in self.networks:
                        scale_model['networks'][net_type] = self.networks[net_type]
                
                # 세포 수준 매개변수 설정
                scale_model['parameters'] = {
                    'cell_division_rate': 0.01,
                    'apoptosis_rate': 0.005,
                    'gene_expression_rates': {'default': 0.2}
                }
                
                multiscale_model['scales']['cellular'] = scale_model
            
            elif scale == 'tissue':
                # 조직 수준 모델 (세포 네트워크, 미세환경 등)
                scale_model = {
                    'structure': 'grid_2d',  # 2D 격자
                    'dimensions': [10, 10],  # 10x10 격자
                    'parameters': {
                        'cell_density': 0.8,
                        'ecm_density': 0.2,
                        'diffusion_constants': {
                            'oxygen': 0.1,
                            'glucose': 0.05,
                            'waste': 0.02
                        }
                    }
                }
                
                multiscale_model['scales']['tissue'] = scale_model
            
            elif scale == 'organism':
                # 유기체 수준 모델 (조직 상호작용, 약물 동태학 등)
                scale_model = {
                    'compartments': ['blood', 'liver', 'kidney', 'tumor'],
                    'parameters': {
                        'blood_flow_rates': {
                            'liver': 0.3,
                            'kidney': 0.2,
                            'tumor': 0.05
                        },
                        'drug_clearance_rates': {
                            'liver': 0.01,
                            'kidney': 0.02
                        }
                    }
                }
                
                multiscale_model['scales']['organism'] = scale_model
        
        # 스케일 간 연결 설정
        if integration_method == 'vertical':
            # 수직 통합 (상위 스케일이 하위 스케일에 의존)
            if 'molecular' in model_scales and 'cellular' in model_scales:
                multiscale_model['scale_connections'].append({
                    'from': 'molecular',
                    'to': 'cellular',
                    'type': 'bottom_up',
                    'mapping': {
                        'metabolites_to_gene_expression': 'function_1',
                        'protein_activity_to_signaling': 'function_2'
                    }
                })
            
            if 'cellular' in model_scales and 'tissue' in model_scales:
                multiscale_model['scale_connections'].append({
                    'from': 'cellular',
                    'to': 'tissue',
                    'type': 'bottom_up',
                    'mapping': {
                        'gene_expression_to_cell_behavior': 'function_3',
                        'signaling_to_cell_movement': 'function_4'
                    }
                })
            
            if 'tissue' in model_scales and 'organism' in model_scales:
                multiscale_model['scale_connections'].append({
                    'from': 'tissue',
                    'to': 'organism',
                    'type': 'bottom_up',
                    'mapping': {
                        'tissue_function_to_organ_function': 'function_5',
                        'local_concentration_to_systemic': 'function_6'
                    }
                })
        
        elif integration_method == 'horizontal':
            # 수평 통합 (스케일 간 대등한 상호작용)
            for i, scale1 in enumerate(model_scales):
                for scale2 in model_scales[i+1:]:
                    multiscale_model['scale_connections'].append({
                        'from': scale1,
                        'to': scale2,
                        'type': 'bidirectional',
                        'mapping': {
                            f"{scale1}_to_{scale2}": f"function_{scale1}_{scale2}",
                            f"{scale2}_to_{scale1}": f"function_{scale2}_{scale1}"
                        }
                    })
        
        elif integration_method == 'nested':
            # 중첩 통합 (상위 스케일 내에 하위 스케일 포함)
            for i in range(len(model_scales) - 1):
                multiscale_model['scale_connections'].append({
                    'from': model_scales[i],
                    'to': model_scales[i+1],
                    'type': 'nested',
                    'embedding': {
                        'scale_ratio': 0.1,  # 스케일 비율
                        'update_frequency': 10  # 업데이트 빈도
                    }
                })
        
        # 다중 스케일 모델 저장
        self.multiscale['model'] = multiscale_model
        
        return multiscale_model
    
    def run_multiscale_simulation(self, steps=10, save_interval=1):
        """
        다중 스케일 모델 시뮬레이션
        
        Parameters:
        -----------
        steps : int, optional
            시뮬레이션 단계 수
        save_interval : int, optional
            결과 저장 간격
        """
        if 'model' not in self.multiscale:
            raise ValueError("먼저 다중 스케일 모델을 생성하세요.")
        
        model = self.multiscale['model']
        scales = model['scales']
        integration_method = model['integration_method']
        
        # 시뮬레이션 초기화
        simulation_results = {scale: [] for scale in scales}
        
        # 간단한 시뮬레이션 구현 (예시)
        for step in range(steps + 1):
            print(f"시뮬레이션 단계 {step}/{steps}")
            
            # 각 스케일별 업데이트
            for scale, scale_model in scales.items():
                if scale == 'molecular':
                    # 분자 수준 시뮬레이션 (간단한 예시)
                    if step == 0:
                        # 초기 상태 (무작위)
                        state = {
                            'metabolite_concentrations': {f"m{i}": np.random.rand() for i in range(10)},
                            'protein_activities': {f"p{i}": np.random.rand() for i in range(10)}
                        }
                    else:
                        # 이전 상태 업데이트
                        prev_state = simulation_results[scale][-1]
                        
                        # 간단한 업데이트 규칙 (확산)
                        diffusion_rate = scale_model['parameters']['diffusion_rate']
                        
                        state = {
                            'metabolite_concentrations': {},
                            'protein_activities': {}
                        }
                        
                        for m, conc in prev_state['metabolite_concentrations'].items():
                            # 무작위 변동 추가
                            new_conc = conc + diffusion_rate * (np.random.rand() - 0.5)
                            state['metabolite_concentrations'][m] = max(0, min(1, new_conc))
                        
                        for p, activity in prev_state['protein_activities'].items():
                            # 무작위 변동 추가
                            new_activity = activity + diffusion_rate * (np.random.rand() - 0.5)
                            state['protein_activities'][p] = max(0, min(1, new_activity))
                    
                    # 결과 저장 (저장 간격에 따라)
                    if step % save_interval == 0:
                        simulation_results[scale].append(state)
                
                elif scale == 'cellular':
                    # 세포 수준 시뮬레이션 (간단한 예시)
                    if step == 0:
                        # 초기 상태 (무작위)
                        state = {
                            'gene_expression': {f"g{i}": np.random.rand() for i in range(10)},
                            'signaling_activities': {f"s{i}": np.random.rand() for i in range(5)},
                            'cell_count': 100
                        }
                    else:
                        # 이전 상태 업데이트
                        prev_state = simulation_results[scale][-1]
                        
                        # 분자 수준의 정보 통합 (있을 경우)
                        if integration_method == 'vertical' and 'molecular' in scales and step > 0:
                            molecular_state = simulation_results['molecular'][-1]
                            
                            # 분자 수준 상태가 세포 수준에 미치는 영향 (단순한 예시)
                            avg_metabolite = np.mean(list(molecular_state['metabolite_concentrations'].values()))
                            avg_protein = np.mean(list(molecular_state['protein_activities'].values()))
                        else:
                            avg_metabolite = 0.5
                            avg_protein = 0.5
                        
                        # 간단한 업데이트 규칙
                        state = {
                            'gene_expression': {},
                            'signaling_activities': {},
                            'cell_count': 0
                        }
                        
                        for g, expr in prev_state['gene_expression'].items():
                            # 대사물질 농도에 영향을 받는 유전자 발현
                            new_expr = expr * (0.8 + 0.4 * avg_metabolite)
                            state['gene_expression'][g] = max(0, min(1, new_expr))
                        
                        for s, activity in prev_state['signaling_activities'].items():
                            # 단백질 활성도에 영향을 받는 신호 활성
                            new_activity = activity * (0.8 + 0.4 * avg_protein)
                            state['signaling_activities'][s] = max(0, min(1, new_activity))
                        
                        # 세포 수 업데이트
                        cell_division_rate = scale_model['parameters']['cell_division_rate']
                        apoptosis_rate = scale_model['parameters']['apoptosis_rate']
                        
                        # 세포 증식 및 사멸
                        growth_rate = cell_division_rate * avg_metabolite - apoptosis_rate
                        state['cell_count'] = max(0, prev_state['cell_count'] * (1 + growth_rate))
                    
                    # 결과 저장 (저장 간격에 따라)
                    if step % save_interval == 0:
                        simulation_results[scale].append(state)
                
                elif scale == 'tissue':
                    # 조직 수준 시뮬레이션 (간단한 예시)
                    if step == 0:
                        # 초기 상태 (격자 기반)
                        dimensions = scale_model['dimensions']
                        grid = np.random.rand(dimensions[0], dimensions[1])
                        state = {
                            'cell_grid': grid,
                            'oxygen_concentration': np.ones(dimensions) * 0.8,
                            'glucose_concentration': np.ones(dimensions) * 0.7
                        }
                    else:
                        # 이전 상태 업데이트
                        prev_state = simulation_results[scale][-1]
                        dimensions = scale_model['dimensions']
                        
                        # 세포 수준의 정보 통합 (있을 경우)
                        if integration_method == 'vertical' and 'cellular' in scales and step > 0:
                            cellular_state = simulation_results['cellular'][-1]
                            avg_gene_expr = np.mean(list(cellular_state['gene_expression'].values()))
                            cell_count = cellular_state['cell_count']
                        else:
                            avg_gene_expr = 0.5
                            cell_count = 100
                        
                        # 간단한 업데이트 규칙 (확산)
                        diffusion_constants = scale_model['parameters']['diffusion_constants']
                        
                        # 이전 상태 복사
                        cell_grid = prev_state['cell_grid'].copy()
                        oxygen = prev_state['oxygen_concentration'].copy()
                        glucose = prev_state['glucose_concentration'].copy()
                        
                        # 세포 성장 및 확산 (단순화된 계산)
                        for i in range(dimensions[0]):
                            for j in range(dimensions[1]):
                                # 산소 및 포도당 확산
                                oxygen_diff = diffusion_constants['oxygen'] * (np.random.rand() - 0.5)
                                glucose_diff = diffusion_constants['glucose'] * (np.random.rand() - 0.5)
                                
                                oxygen[i, j] += oxygen_diff
                                glucose[i, j] += glucose_diff
                                
                                # 경계 조건
                                oxygen[i, j] = max(0, min(1, oxygen[i, j]))
                                glucose[i, j] = max(0, min(1, glucose[i, j]))
                                
                                # 세포 성장률 (산소 및 포도당 농도에 비례)
                                growth_factor = (oxygen[i, j] + glucose[i, j]) / 2
                                cell_grid[i, j] += 0.1 * growth_factor * (avg_gene_expr - 0.5)
                                
                                # 경계 조건
                                cell_grid[i, j] = max(0, min(1, cell_grid[i, j]))
                        
                        # 세포 수와 격자 밀도 일치
                        total_density = np.sum(cell_grid) / (dimensions[0] * dimensions[1])
                        scale_factor = cell_count / 100 / total_density
                        cell_grid *= scale_factor
                        
                        state = {
                            'cell_grid': cell_grid,
                            'oxygen_concentration': oxygen,
                            'glucose_concentration': glucose
                        }
                    
                    # 결과 저장 (저장 간격에 따라)
                    if step % save_interval == 0:
                        simulation_results[scale].append(state)
                
                elif scale == 'organism':
                    # 유기체 수준 시뮬레이션 (간단한 예시)
                    if step == 0:
                        # 초기 상태
                        state = {
                            'compartment_concentrations': {
                                'blood': 1.0,  # 약물 초기 농도
                                'liver': 0.0,
                                'kidney': 0.0,
                                'tumor': 0.0
                            },
                            'physiological_parameters': {
                                'blood_pressure': 120,
                                'heart_rate': 70,
                                'body_temperature': 37.0
                            }
                        }
                    else:
                        # 이전 상태 업데이트
                        prev_state = simulation_results[scale][-1]
                        
                        # 조직 수준의 정보 통합 (있을 경우)
                        if integration_method == 'vertical' and 'tissue' in scales and step > 0:
                            tissue_state = simulation_results['tissue'][-1]
                            avg_cell_density = np.mean(tissue_state['cell_grid'])
                        else:
                            avg_cell_density = 0.5
                        
                        # 약물 분포 업데이트 (간단한 구획 모델)
                        concentrations = prev_state['compartment_concentrations'].copy()
                        flow_rates = scale_model['parameters']['blood_flow_rates']
                        clearance_rates = scale_model['parameters']['drug_clearance_rates']
                        
                        # 혈액에서 각 조직으로의 이동
                        blood_conc = concentrations['blood']
                        
                        # 각 구획으로의 분포
                        for comp, flow in flow_rates.items():
                            transfer = blood_conc * flow * 0.1
                            concentrations['blood'] -= transfer
                            concentrations[comp] += transfer
                        
                        # 약물 제거
                        for comp, clear in clearance_rates.items():
                            removal = concentrations[comp] * clear
                            concentrations[comp] -= removal
                        
                        # 생리학적 매개변수 업데이트 (약물 농도에 따라)
                        phys_params = prev_state['physiological_parameters'].copy()
                        
                        # 약물이 혈압에 미치는 영향 (간단한 예시)
                        if blood_conc > 0.5:
                            phys_params['blood_pressure'] -= 5
                        
                        # 종양 크기가 체온에 미치는 영향 (간단한 예시)
                        tumor_effect = concentrations['tumor'] * avg_cell_density
                        phys_params['body_temperature'] += 0.1 * tumor_effect
                        
                        state = {
                            'compartment_concentrations': concentrations,
                            'physiological_parameters': phys_params
                        }
                    
                    # 결과 저장 (저장 간격에 따라)
                    if step % save_interval == 0:
                        simulation_results[scale].append(state)
        
        # 시뮬레이션 결과 저장
        self.multiscale['simulation_results'] = simulation_results
        self.multiscale['simulation_steps'] = steps
        self.multiscale['save_interval'] = save_interval
        
        return simulation_results
    
    def visualize_multiscale_simulation(self, scale=None):
        """
        다중 스케일 시뮬레이션 결과 시각화
        
        Parameters:
        -----------
        scale : str, optional
            시각화할 특정 스케일 (None이면 모든 스케일)
        """
        if 'simulation_results' not in self.multiscale:
            raise ValueError("먼저 다중 스케일 시뮬레이션을 실행하세요.")
        
        results = self.multiscale['simulation_results']
        steps = self.multiscale['simulation_steps']
        save_interval = self.multiscale['save_interval']
        
        # 시각화할 스케일 선택
        if scale is not None:
            if scale not in results:
                raise ValueError(f"'{scale}' 스케일을 찾을 수 없습니다.")
            
            scales_to_visualize = [scale]
        else:
            scales_to_visualize = list(results.keys())
        
        # 각 스케일별 시각화
        for scale in scales_to_visualize:
            scale_results = results[scale]
            time_points = list(range(0, steps + 1, save_interval))
            
            if scale == 'molecular':
                # 분자 수준 시각화
                
                # 1. 대사물질 농도 변화
                metabolites = list(scale_results[0]['metabolite_concentrations'].keys())
                
                plt.figure(figsize=(12, 8))
                for m in metabolites:
                    concentrations = [result['metabolite_concentrations'][m] for result in scale_results]
                    plt.plot(time_points, concentrations, marker='o', label=m)
                
                plt.xlabel('Simulation Step')
                plt.ylabel('Concentration')
                plt.title('Metabolite Concentrations Over Time')
                plt.legend()
                plt.grid(True, linestyle='--', alpha=0.7)
                
                # 결과 저장
                plt.savefig(os.path.join(self.results_dir, "molecular_metabolites.png"), dpi=300, bbox_inches='tight')
                plt.close()
                
                # 2. 단백질 활성도 변화
                proteins = list(scale_results[0]['protein_activities'].keys())
                
                plt.figure(figsize=(12, 8))
                for p in proteins:
                    activities = [result['protein_activities'][p] for result in scale_results]
                    plt.plot(time_points, activities, marker='o', label=p)
                
                plt.xlabel('Simulation Step')
                plt.ylabel('Activity')
                plt.title('Protein Activities Over Time')
                plt.legend()
                plt.grid(True, linestyle='--', alpha=0.7)
                
                # 결과 저장
                plt.savefig(os.path.join(self.results_dir, "molecular_proteins.png"), dpi=300, bbox_inches='tight')
                plt.close()
            
            elif scale == 'cellular':
                # 세포 수준 시각화
                
                # 1. 유전자 발현 변화
                genes = list(scale_results[0]['gene_expression'].keys())
                
                plt.figure(figsize=(12, 8))
                for g in genes:
                    expressions = [result['gene_expression'][g] for result in scale_results]
                    plt.plot(time_points, expressions, marker='o', label=g)
                
                plt.xlabel('Simulation Step')
                plt.ylabel('Expression Level')
                plt.title('Gene Expression Over Time')
                plt.legend()
                plt.grid(True, linestyle='--', alpha=0.7)
                
                # 결과 저장
                plt.savefig(os.path.join(self.results_dir, "cellular_gene_expression.png"), dpi=300, bbox_inches='tight')
                plt.close()
                
                # 2. 세포 수 변화
                cell_counts = [result['cell_count'] for result in scale_results]
                
                plt.figure(figsize=(10, 6))
                plt.plot(time_points, cell_counts, marker='o', linewidth=2)
                plt.xlabel('Simulation Step')
                plt.ylabel('Cell Count')
                plt.title('Cell Population Over Time')
                plt.grid(True, linestyle='--', alpha=0.7)
                
                # 결과 저장
                plt.savefig(os.path.join(self.results_dir, "cellular_cell_count.png"), dpi=300, bbox_inches='tight')
                plt.close()
            
            elif scale == 'tissue':
                # 조직 수준 시각화
                
                # 1. 세포 분포 히트맵 (선택된 시간 포인트)
                time_indices = [0, len(scale_results)//2, -1]  # 처음, 중간, 마지막
                
                fig, axes = plt.subplots(1, len(time_indices), figsize=(15, 5))
                for i, t_idx in enumerate(time_indices):
                    t = time_points[t_idx]
                    cell_grid = scale_results[t_idx]['cell_grid']
                    
                    im = axes[i].imshow(cell_grid, cmap='viridis', vmin=0, vmax=1)
                    axes[i].set_title(f'Cell Distribution at t={t}')
                    axes[i].axis('off')
                
                # 컬러바 추가
                cbar = fig.colorbar(im, ax=axes, orientation='horizontal', pad=0.05)
                cbar.set_label('Cell Density')
                
                plt.tight_layout()
                
                # 결과 저장
                plt.savefig(os.path.join(self.results_dir, "tissue_cell_distribution.png"), dpi=300, bbox_inches='tight')
                plt.close()
                
                # 2. 산소 농도 시간 변화 (평균)
                oxygen_means = [np.mean(result['oxygen_concentration']) for result in scale_results]
                glucose_means = [np.mean(result['glucose_concentration']) for result in scale_results]
                
                plt.figure(figsize=(10, 6))
                plt.plot(time_points, oxygen_means, marker='o', linewidth=2, label='Oxygen')
                plt.plot(time_points, glucose_means, marker='s', linewidth=2, label='Glucose')
                plt.xlabel('Simulation Step')
                plt.ylabel('Average Concentration')
                plt.title('Nutrient Concentration Over Time')
                plt.legend()
                plt.grid(True, linestyle='--', alpha=0.7)
                
                # 결과 저장
                plt.savefig(os.path.join(self.results_dir, "tissue_nutrients.png"), dpi=300, bbox_inches='tight')
                plt.close()
            
            elif scale == 'organism':
                # 유기체 수준 시각화
                
                # 1. 구획별 약물 농도 변화
               # 1. 구획별 약물 농도 변화
                compartments = list(scale_results[0]['compartment_concentrations'].keys())
                
                plt.figure(figsize=(12, 8))
                for comp in compartments:
                    concentrations = [result['compartment_concentrations'][comp] for result in scale_results]
                    plt.plot(time_points, concentrations, marker='o', linewidth=2, label=comp.capitalize())
                
                plt.xlabel('Simulation Step')
                plt.ylabel('Drug Concentration')
                plt.title('Drug Distribution Across Compartments')
                plt.legend()
                plt.grid(True, linestyle='--', alpha=0.7)
                
                # 결과 저장
                plt.savefig(os.path.join(self.results_dir, "organism_drug_distribution.png"), dpi=300, bbox_inches='tight')
                plt.close()
                
                # 2. 생리학적 매개변수 변화
                params = list(scale_results[0]['physiological_parameters'].keys())
                
                fig, axes = plt.subplots(len(params), 1, figsize=(12, 4*len(params)), sharex=True)
                
                for i, param in enumerate(params):
                    values = [result['physiological_parameters'][param] for result in scale_results]
                    axes[i].plot(time_points, values, marker='o', linewidth=2)
                    axes[i].set_ylabel(param.replace('_', ' ').title())
                    axes[i].grid(True, linestyle='--', alpha=0.7)
                
                axes[-1].set_xlabel('Simulation Step')
                fig.suptitle('Physiological Parameters Over Time')
                plt.tight_layout()
                plt.subplots_adjust(top=0.95)
                
                # 결과 저장
                plt.savefig(os.path.join(self.results_dir, "organism_physiological_parameters.png"), dpi=300, bbox_inches='tight')
                plt.close()
        
        # 다중 스케일 통합 시각화 (스케일 간 관계)
        if len(scales_to_visualize) > 1 and self.multiscale['model']['integration_method'] == 'vertical':
            # 수직 통합 모델에서 스케일 간 관계 시각화
            
            # 스케일 간 중요 지표 선택
            key_metrics = {}
            
            if 'molecular' in scales_to_visualize:
                # 분자 수준 주요 지표 (대사물질 및 단백질 평균)
                metabolite_means = []
                protein_means = []
                
                for result in results['molecular']:
                    metabolite_means.append(np.mean(list(result['metabolite_concentrations'].values())))
                    protein_means.append(np.mean(list(result['protein_activities'].values())))
                
                key_metrics['molecular'] = {
                    'metabolite_mean': metabolite_means,
                    'protein_mean': protein_means
                }
            
            if 'cellular' in scales_to_visualize:
                # 세포 수준 주요 지표 (유전자 발현 및 세포 수)
                gene_means = []
                cell_counts = []
                
                for result in results['cellular']:
                    gene_means.append(np.mean(list(result['gene_expression'].values())))
                    cell_counts.append(result['cell_count'])
                
                key_metrics['cellular'] = {
                    'gene_expression_mean': gene_means,
                    'cell_count': cell_counts
                }
            
            if 'tissue' in scales_to_visualize:
                # 조직 수준 주요 지표 (세포 밀도 및 산소 농도)
                cell_density_means = []
                oxygen_means = []
                
                for result in results['tissue']:
                    cell_density_means.append(np.mean(result['cell_grid']))
                    oxygen_means.append(np.mean(result['oxygen_concentration']))
                
                key_metrics['tissue'] = {
                    'cell_density_mean': cell_density_means,
                    'oxygen_mean': oxygen_means
                }
            
            if 'organism' in scales_to_visualize:
                # 유기체 수준 주요 지표 (종양 구획 약물 농도 및 체온)
                tumor_conc = []
                body_temp = []
                
                for result in results['organism']:
                    tumor_conc.append(result['compartment_concentrations'].get('tumor', 0))
                    body_temp.append(result['physiological_parameters'].get('body_temperature', 37))
                
                key_metrics['organism'] = {
                    'tumor_drug_concentration': tumor_conc,
                    'body_temperature': body_temp
                }
            
            # 스케일 간 관계 시각화
            for i, scale1 in enumerate(scales_to_visualize[:-1]):
                for scale2 in scales_to_visualize[i+1:]:
                    # 두 스케일 간의 모든 지표 쌍에 대한 상관관계 시각화
                    metrics1 = key_metrics.get(scale1, {})
                    metrics2 = key_metrics.get(scale2, {})
                    
                    if not metrics1 or not metrics2:
                        continue
                    
                    # 그리드 크기 계산
                    num_metrics1 = len(metrics1)
                    num_metrics2 = len(metrics2)
                    
                    fig, axes = plt.subplots(num_metrics1, num_metrics2, figsize=(4*num_metrics2, 4*num_metrics1))
                    
                    # 단일 지표인 경우 2D 배열로 변환
                    if num_metrics1 == 1 and num_metrics2 == 1:
                        axes = np.array([[axes]])
                    elif num_metrics1 == 1:
                        axes = axes.reshape(1, -1)
                    elif num_metrics2 == 1:
                        axes = axes.reshape(-1, 1)
                    
                    for row, (metric1_name, metric1_values) in enumerate(metrics1.items()):
                        for col, (metric2_name, metric2_values) in enumerate(metrics2.items()):
                            ax = axes[row, col]
                            
                            # 산점도 그리기
                            ax.scatter(metric1_values, metric2_values, alpha=0.7)
                            
                            # 추세선 추가
                            try:
                                z = np.polyfit(metric1_values, metric2_values, 1)
                                p = np.poly1d(z)
                                ax.plot(metric1_values, p(metric1_values), "r--", alpha=0.7)
                                
                                # 상관 계수 계산
                                correlation = np.corrcoef(metric1_values, metric2_values)[0, 1]
                                ax.text(0.05, 0.95, f"r = {correlation:.2f}", transform=ax.transAxes,
                                        verticalalignment='top', bbox=dict(boxstyle='round', alpha=0.5))
                            except:
                                pass
                            
                            # 축 레이블 설정
                            ax.set_xlabel(metric1_name.replace('_', ' ').title())
                            ax.set_ylabel(metric2_name.replace('_', ' ').title())
                            ax.grid(True, linestyle='--', alpha=0.7)
                    
                    fig.suptitle(f'Relationships Between {scale1.capitalize()} and {scale2.capitalize()} Scales')
                    plt.tight_layout()
                    plt.subplots_adjust(top=0.95)
                    
                    # 결과 저장
                    plt.savefig(os.path.join(self.results_dir, f"{scale1}_{scale2}_relationships.png"), dpi=300, bbox_inches='tight')
                    plt.close()
# OMNIDRUG-AI Platform

AI 기반 신약개발 통합 플랫폼 - 프로젝트 단위 파이프라인 관리 시스템

## 📋 목차
- [개요](#개요)
- [주요 기능](#주요-기능)
- [기술 스택](#기술-스택)
- [설치 방법](#설치-방법)
- [실행 방법](#실행-방법)
- [API 문서](#api-문서)
- [프로젝트 구조](#프로젝트-구조)
- [개발 현황](#개발-현황)

## 🔬 개요

OMNIDRUG-AI는 인공지능과 컴퓨터 과학을 활용한 신약개발 플랫폼입니다. 가상 스크리닝, 분자 도킹, ADMET 예측 등의 기능을 통합하여 신약개발 프로세스를 가속화합니다.

## ✨ 주요 기능

- **프로젝트 관리**: 타겟 단백질 기반 신약개발 프로젝트 관리
- **가상 스크리닝**: UniDock을 활용한 분자 도킹
- **분자 데이터 관리**: SMILES, InChI 기반 화합물 데이터베이스
- **배치 파이프라인**: 스크리닝-합성-분석 워크플로우
- **실시간 협업**: 다중 사용자 프로젝트 관리

## 🛠️ 기술 스택

### Backend
- Python 3.11+
- Django 5.0.1
- Django REST Framework
- Celery + Redis (비동기 작업)
- PostgreSQL / SQLite3
- RDKit (화학정보학)

### Frontend (개발 중)
- Vue.js 3
- Vite
- Element Plus
- Pinia (상태관리)
- ECharts (시각화)

## 📦 설치 방법

### 1. 저장소 클론
```bash
git clone https://github.com/SONGINHYEOK/CODE.git
cd CODE/OMNIDRUG-AI
```

### 2. Python 환경 설정 (Conda 권장)
```bash
# Conda 환경 생성
conda create -n omnidrug python=3.11
conda activate omnidrug

# RDKit 설치 (conda-forge 채널)
conda install -c conda-forge rdkit
```

### 3. Python 패키지 설치
```bash
cd backend
pip install -r requirements.txt
```

### 4. 환경 변수 설정
```bash
# .env 파일 생성
cp .env.example .env

# .env 파일 편집
# SECRET_KEY, DATABASE 설정 등 수정
```

### 5. 데이터베이스 초기화
```bash
# 마이그레이션 생성
python manage.py makemigrations core
python manage.py makemigrations projects
python manage.py makemigrations molecules
python manage.py makemigrations screening

# 마이그레이션 적용
python manage.py migrate

# 슈퍼유저 생성
python manage.py createsuperuser
```

## 🚀 실행 방법

### Backend 서버 실행
```bash
# 1. Conda 환경 활성화
conda activate omnidrug

# 2. backend 디렉토리로 이동
cd backend

# 3. Django 개발 서버 실행
python manage.py runserver

# 서버가 http://127.0.0.1:8000 에서 실행됩니다
```

### Celery Worker 실행 (별도 터미널)
```bash
# Redis 먼저 실행 (Mac)
redis-server

# Celery worker 실행
cd backend
celery -A omnidrug worker -l info

# Flower (Celery 모니터링) - 선택사항
celery -A omnidrug flower
```

### Frontend 실행 (개발 중)
```bash
cd frontend
npm install
npm run dev

# http://localhost:5173 에서 실행
```

## 📚 API 문서

서버 실행 후 다음 URL에서 API 문서를 확인할 수 있습니다:

- **Swagger UI**: http://127.0.0.1:8000/api/docs/
- **Django Admin**: http://127.0.0.1:8000/admin/

### 주요 API 엔드포인트

#### 인증
- `POST /api/auth/login/` - JWT 토큰 로그인
- `POST /api/auth/register/` - 회원가입
- `GET /api/auth/profile/` - 프로필 조회

#### 프로젝트 관리
- `GET /api/projects/` - 프로젝트 목록
- `POST /api/projects/` - 프로젝트 생성
- `GET /api/projects/{id}/` - 프로젝트 상세
- `PUT /api/projects/{id}/` - 프로젝트 수정
- `DELETE /api/projects/{id}/` - 프로젝트 삭제

#### 타겟 관리
- `GET /api/targets/` - 타겟 목록
- `POST /api/targets/` - 타겟 추가

#### 배치 관리
- `GET /api/batches/` - 배치 목록
- `POST /api/batches/` - 배치 생성

## 📁 프로젝트 구조

```
OMNIDRUG-AI/
├── backend/
│   ├── omnidrug/           # Django 프로젝트 설정
│   │   ├── settings.py     # 설정 파일
│   │   ├── urls.py         # URL 라우팅
│   │   ├── celery.py       # Celery 설정
│   │   └── wsgi.py         # WSGI 설정
│   │
│   ├── apps/               # Django 앱
│   │   ├── core/          # 사용자 인증
│   │   ├── projects/      # 프로젝트 관리
│   │   ├── molecules/     # 분자 데이터
│   │   └── screening/     # 가상 스크리닝
│   │
│   ├── requirements.txt    # Python 패키지
│   ├── manage.py          # Django 관리 스크립트
│   └── db.sqlite3         # 개발용 데이터베이스
│
├── frontend/              # Vue.js 프론트엔드
│   ├── src/
│   ├── public/
│   └── package.json
│
├── docker/                # Docker 설정
├── scripts/               # 유틸리티 스크립트
├── .env.example          # 환경변수 예제
├── .gitignore
└── README.md
```

## 🧪 테스트 실행

```bash
# Django 테스트
python manage.py test

# 특정 앱 테스트
python manage.py test apps.projects
```

## 🐳 Docker 실행 (선택사항)

```bash
# Docker Compose로 전체 스택 실행
docker-compose up -d

# 로그 확인
docker-compose logs -f
```

## 📊 개발 현황

### 완료된 기능 ✅
- Django REST API 구축
- JWT 인증 시스템
- 프로젝트/타겟/배치 CRUD
- Swagger API 문서화
- Celery 비동기 작업 설정

### 진행 중 🚧
- Vue.js 프론트엔드 개발
- UniDock 통합
- 분자 시각화 기능

### 예정된 기능 📝
- ADMET 예측 모델
- 클러스터링 분석
- SAR 분석
- 실시간 알림 시스템

## 🤝 기여 방법

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📝 라이선스

이 프로젝트는 MIT 라이선스 하에 있습니다.

## 👥 개발팀

- **송인혁** - 프로젝트 리드

## 📞 문의

프로젝트 관련 문의사항은 Issues 탭을 이용해주세요.

---

## ⚡ Quick Start

빠른 시작을 위한 명령어 모음:

```bash
# 1. 클론 및 이동
git clone https://github.com/SONGINHYEOK/CODE.git
cd CODE/OMNIDRUG-AI/backend

# 2. Python 환경 (이미 설치되어 있다면 활성화만)
conda activate omnidrug

# 3. 패키지 설치 (최초 1회)
pip install -r requirements.txt

# 4. DB 초기화 (최초 1회)
python manage.py migrate
python manage.py createsuperuser

# 5. 서버 실행
python manage.py runserver

# 6. 브라우저에서 접속
# http://127.0.0.1:8000/api/docs/
```

## 🔧 문제 해결

### 자주 발생하는 문제들

#### 1. ModuleNotFoundError
```bash
# 해결: 패키지 재설치
pip install -r requirements.txt
```

#### 2. Database 에러
```bash
# 해결: 데이터베이스 초기화
rm db.sqlite3
python manage.py migrate
```

#### 3. Port 8000 already in use
```bash
# 해결: 다른 포트 사용
python manage.py runserver 8001
```

---

*Last Updated: 2024-12-20*
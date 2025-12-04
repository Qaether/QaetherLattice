
# 루프‑위상 스캐너 (Loop-Phase Scanner, FCC)

이 템플릿은 **포논 탑올로지 DB**에서 내보낸(또는 Phonopy/DFPT에서 얻은) 고유벡터/주파수 자료를
불러와 **FCC 격자의 삼각/사각 최소 루프**에서 밴드별 **Wilson loop 위상**을 수치적으로 적분하고,
각 위상이 `π/6`의 정수배인지 판정합니다.

## 빠른 시작

```bash
# 가상환경 권장
pip install numpy pyyaml

# 설정 파일 수정 (example_config.yaml)
#  - input_npz 를 실제 데이터 경로로 교체
#  - bands, tolerance 조정 가능

python loop_phase_scanner.py --config example_config.yaml
```

출력: `scan_results.json` 에 각 루프/밴드의 위상 및 `π/6` 정수배 여부가 기록됩니다.

## 입력 데이터 포맷 (.npz)

키 | 형태 | 설명
---|---|---
`kpts` | `(Nk,3)` | 역격자 k-점 (2π/a 단위 또는 crystal, 자유)
`evals`| `(Nk,Nb)` | 주파수(THz) 또는 ω² (옵션)
`evecs`| `(Nk,Nb,Ndof)` | 복소 고유벡터 (질량가중 정규화 권장)
`grid_shape`| `(3,)` | 정수 튜플 (nx,ny,nz), 직교 grid 가정

> **주의:** 게이지 연속성(위상)가 깨지면 Wilson loop가 잡음이 커집니다.
> 가능하면 **Wannier‑interpolation phonon** 또는 **연속 게이지 고정**을 권장합니다.

## 데이터 획득 팁

- **Topological Phonon Database (TPDB)**: 이상적 후보군 선별 후, 표면/대역 자료를 CSV/JSON으로 내보낸 다음
  본 스캐너에 맞춰 변환하세요. (TPDB 소개 및 사용례: Science 384, eadf8458 (2024)).
- **NIST JARVIS‑DFT**: DFPT 기반 포논 및 원시 입/출력 제공. (phonon DOS/DFPT raw I/O 데이터셋 참조).

## 매개변수

- `unit` : 검증 단위 (기본 `π/6`)
- `tolerance` : 허용 오차 (라디안)
- `bands` : 검사할 밴드 범위 또는 리스트
- `max_loops` : 과도한 계산 방지용 샘플링 한도

## 방법 개요

Wilson loop 위상은 이웃 k점 고유벡터의 병렬 수송 내적의 복소 위상을 합산하여 계산합니다.
FCC의 짧은 **삼각(111)** 및 **사각(100/110)** 루프에서 `3Δφ`와 `4Δφ`가 각각 `2πℤ`를 강제하므로,
공통 최소 단위는 `gcd(2π/3, π/2)=π/6` 입니다. (수학적 배경은 본 프로젝트 문서의 양자화 노트 참조).

## 결과 해석

- `phase` : 루프의 총 위상 (주값 브랜치, `(-π,π]`)
- `nearest_m`: 가장 가까운 정수 `m` (즉, `m * unit ≈ phase`)
- `delta` : 편차
- `ok_pi_over_6` : `|delta| ≤ tolerance` 이면 `True`

## 라이선스

연구/프로토타입 용도. 상용 사용 전 검증 필요.

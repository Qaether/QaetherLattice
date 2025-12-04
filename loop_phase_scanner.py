
"""
루프‑위상 스캐너 (Loop-Phase Scanner for Phonon Topology, FCC 전용)
Author: Qaether v1.3
Python >= 3.9

개요
----
- 입력: k-격자 위 phonon 고유벡터(정규직교, 질량가중 기준)와 고유값(주파수), k-점 좌표.
- 처리: FCC 실격자 ↔ BCC 역격자에서 (111) 삼각 루프, (100)/(110) 사각 루프를 생성하여
        각 밴드에 대한 Wilson loop 위상 (Berry phase) 를 수치적으로 적분.
- 출력: 각 루프/밴드의 위상합이 π/6 의 정수배인지(허용 오차 내) 판정, 히트맵/표 CSV.

데이터 형태(권장)
----------------
npz 파일에 다음 키를 포함:
- kpts : (Nk,3) reciprocal coords in Cartesian 2π/a units (또는 crystal coords, config에서 지정)
- evals: (Nk, Nb) phonon 주파수(THz) 또는 ω^2 (단위 자유)
- evecs: (Nk, Nb, Ndof) 복소 고유벡터, Ndof=3Natoms (질량가중 정규화 권장)
- grid_shape: (nx,ny,nz) 정수 튜플 (직교 grid로 가정)
주의: 위상 적분의 gauge-smoothness 확보를 위해 각 k점에서 밴드별 위상 gauge가 연속적이어야 함.
      실데이터는 Wannier 기반(phonon-wannier), 또는 DFPT 인터폴레이션 결과를 권장.

사용법
------
python loop_phase_scanner.py --config example_config.yaml
"""
import yaml
import argparse, json, math
from pathlib import Path
import numpy as np

# ------------------------ 유틸 ------------------------

def _complex_overlap(u, v):
    """밴드 스피너(고유벡터) 내적: <u|v>, returns complex scalar."""
    return np.vdot(u, v)  # mass-weighted 정규화가 입력에서 이미 되어있다고 가정

def _link_phase(u, v, clip=1.0):
    """이웃 k점 사이의 링크 위상(게이지-인자 제거용 병렬 수송)"""
    ov = _complex_overlap(u, v)
    # 안정화를 위해 정규화 (수치 오차로 |ov|>1 약간 넘어갈 수 있음)
    x = np.clip(np.real(ov), -clip, clip) + 1j*np.clip(np.imag(ov), -clip, clip)
    return np.angle(x)

def _wilson_loop_phase(evecs_loop):
    """
    evecs_loop: list of eigenvectors for a fixed band at ordered loop vertices [k0,k1,...,kN(=k0)]
    Wilson loop 위상 = 링크 위상의 합 (병렬 수송 gauge)
    """
    total = 0.0
    for a, b in zip(evecs_loop[:-1], evecs_loop[1:]):
        total += _link_phase(a, b)
    # 2π wrap to principal branch
    return ( (total + np.pi) % (2*np.pi) ) - np.pi

def _closest_multiple(x, unit):
    """x가 unit의 정수배에 얼마나 가까운지 반환 (정수 m, 오차 delta)"""
    m = int(round(x / unit))
    return m, x - m*unit

# ------------------------ 루프 구축 (BCC reciprocal) ------------------------

def _grid_index(ix, iy, iz, shape):
    nx, ny, nz = shape
    return ((ix % nx) * ny + (iy % ny)) * nz + (iz % nz)

def _collect_band_vectors(evecs, loop_indices, band):
    """evecs: (Nk,Nb,Ndof), loop_indices: list[int]"""
    return [evecs[idx, band] for idx in loop_indices]

def _make_square_loop(ix, iy, iz, ax, ay, shape):
    """
    정사각 루프: 시작 (ix,iy,iz)에서 x-축 ax, y-축 ay 방향 한 칸씩
    반환: 인덱스 리스트 (마지막에 시작점으로 닫힘)
    """
    path = [
        (ix, iy, iz),
        (ix+ax, iy, iz),
        (ix+ax, iy+ay, iz),
        (ix, iy+ay, iz),
        (ix, iy, iz),
    ]
    return [_grid_index(*p, shape) for p in path]

def _make_triangle_loop(ix, iy, iz, a1, a2, shape):
    """
    삼각 루프: 평면 벡터 a1, a2 (격자 스텝 단위) 사용
    """
    path = [
        (ix, iy, iz),
        (ix+a1[0], iy+a1[1], iz+a1[2]),
        (ix+a2[0], iy+a2[1], iz+a2[2]),
        (ix, iy, iz),
    ]
    return [_grid_index(*p, shape) for p in path]

def generate_loops(shape, mode="fcc_default"):
    """
    shape: (nx,ny,nz)
    mode: "fcc_default" → (100)/(110) 사각 + (111) 삼각 루프 셋
    """
    nx, ny, nz = shape
    loops = []
    # 사각 루프: (x,y) 평면, (x,z), (y,z)
    for ax, ay in [(1,1), (1,0), (0,1)]:
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    if (ax,ay)==(1,1):
                        # (x,y)
                        loops.append(("square_xy", _make_square_loop(ix,iy,iz,1,1,shape)))
                    elif (ax,ay)==(1,0):
                        # (x,z)
                        loops.append(("square_xz", _make_square_loop(ix,iy,iz,1,1, (nx,nz,ny))))
                    else:
                        # (y,z)
                        loops.append(("square_yz", _make_square_loop(ix,iy,iz,1,1, (ny,nz,nx))))
    # 삼각 루프: (111)류 평면 근사 (격자 스텝으로 (1,0,0)->(0,1,0)->(0,0,1))
    tri_steps = [(1,0,0),(0,1,0),(0,0,1)]
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                loops.append(("tri_111", _make_triangle_loop(ix,iy,iz, tri_steps[0], tri_steps[1], shape)))
                loops.append(("tri_111", _make_triangle_loop(ix,iy,iz, tri_steps[1], tri_steps[2], shape)))
                loops.append(("tri_111", _make_triangle_loop(ix,iy,iz, tri_steps[2], tri_steps[0], shape)))
    return loops

# ------------------------ 메인: 위상 스캔 ------------------------

def scan(npz_path, out_json, unit=np.pi/6, tol=0.05, bands=None, max_loops=20000):
    """
    npz_path : 입력 데이터 경로
    out_json : 결과 저장 경로
    unit     : 검증 단위 (기본 π/6)
    tol      : 허용 오차 (라디안)
    bands    : [bmin,bmax) 또는 리스트. None이면 전 밴드
    max_loops: 과도한 계산 방지 (샘플링)
    """
    data = np.load(npz_path, allow_pickle=True)
    evecs = data["evecs"]  # (Nk,Nb,Ndof)
    shape = tuple(data["grid_shape"])
    Nk, Nb, Ndof = evecs.shape
    if bands is None:
        band_list = list(range(Nb))
    elif isinstance(bands, (list, tuple)) and len(bands)==2 and isinstance(bands[0], int):
        band_list = list(range(bands[0], bands[1]))
    else:
        band_list = list(bands)

    loops = generate_loops(shape)
    if len(loops) > max_loops:
        # 균일 샘플링
        step = max(1, len(loops)//max_loops)
        loops = loops[::step]

    results = []
    for lname, lidx in loops:
        for b in band_list:
            vecs = _collect_band_vectors(evecs, lidx, b)
            phi = _wilson_loop_phase(vecs)
            m, delta = _closest_multiple(phi, unit)
            ok = abs(delta) <= tol
            results.append({
                "loop_type": lname,
                "band": int(b),
                "phase": float(phi),
                "nearest_m": int(m),
                "delta": float(delta),
                "ok_pi_over_6": bool(ok)
            })

    # 통계 요약
    hits = sum(1 for r in results if r["ok_pi_over_6"])
    summary = {
        "npz": str(npz_path),
        "Nk": Nk,
        "Nb": Nb,
        "Ndof": Ndof,
        "loops_scanned": len(results),
        "hits_pi_over_6": hits,
        "hit_ratio": hits/max(1,len(results))
    }

    Path(out_json).write_text(json.dumps({"summary": summary, "results": results}, indent=2, ensure_ascii=False))
    return summary

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=str, required=True, help="YAML 설정 파일")
    args = ap.parse_args()
    cfg = yaml.safe_load(Path(args.config).read_text())
    input_data_path= '/Users/francishan/Documents/Github/QaetherLattice/demo_fake_data.npz'
    ouput_data_path= '/Users/francishan/Documents/Github/QaetherLattice/output_json'
    summary = scan(
        npz_path= input_data_path,
        out_json= ouput_data_path,
        unit=cfg.get("unit", math.pi/6),
        tol=cfg.get("tolerance", 0.05),
        bands=cfg.get("bands", None),
        max_loops=cfg.get("max_loops", 20000)
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))

if __name__ == "__main__":
    main()

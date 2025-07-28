import numpy as np
import matplotlib.pyplot as plt

def generate_fcc_lattice(nx: int, ny: int, nz: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Generates unique corner and face-centered atom positions for an FCC lattice.
    FCC 격자에 대한 고유한 꼭짓점 및 면 중심 원자 위치를 생성합니다.

    :param nx: The number of unit cells in the x-dimension.
    :param ny: The number of unit cells in the y-dimension.
    :param nz: The number of unit cells in the z-dimension.
    :return: A tuple containing two numpy arrays: (corners, face_centers).
    """
    # 1. 꼭짓점 원자 생성 (정수 좌표)
    i_c, j_c, k_c = np.mgrid[0:nx + 1, 0:ny + 1, 0:nz + 1]
    corners = np.vstack([i_c.ravel(), j_c.ravel(), k_c.ravel()]).T

    # 2. 면 중심 원자 생성
    # XY 평면의 면 중심
    i_xy, j_xy, k_xy = np.mgrid[0:nx, 0:ny, 0:nz + 1]
    faces_xy = np.vstack([i_xy.ravel() + 0.5, j_xy.ravel() + 0.5, k_xy.ravel()]).T

    # XZ 평면의 면 중심
    i_xz, j_xz, k_xz = np.mgrid[0:nx, 0:ny + 1, 0:nz]
    faces_xz = np.vstack([i_xz.ravel() + 0.5, j_xz.ravel(), k_xz.ravel() + 0.5]).T

    # YZ 평면의 면 중심
    i_yz, j_yz, k_yz = np.mgrid[0:nx + 1, 0:ny, 0:nz]
    faces_yz = np.vstack([i_yz.ravel(), j_yz.ravel() + 0.5, k_yz.ravel() + 0.5]).T

    # 모든 면 중심 원자를 합칩니다.
    face_centers = np.vstack([faces_xy, faces_xz, faces_yz])

    return corners, face_centers

def plot_lattice(corners: np.ndarray, face_centers: np.ndarray, title: str, grid_dims: tuple[int, int, int]):
    """
    Plots the given atomic coordinates in a 3D scatter plot with different colors.
    주어진 원자 좌표를 다른 색상으로 구분하여 3D 산점도로 시각화합니다.
    """
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(projection='3d')

    # 꼭짓점 원자 (파란색)
    ax.scatter(corners[:, 0], corners[:, 1], corners[:, 2], s=60, c='blue', marker='o', alpha=0.8, label='Corners')
    # 면 중심 원자 (빨간색)
    ax.scatter(face_centers[:, 0], face_centers[:, 1], face_centers[:, 2], s=60, c='red', marker='o', alpha=0.8, label='Face Centers')

    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    ax.set_title(title, fontsize=16)

    nx, ny, nz = grid_dims
    ax.set_xlim([-0.2, nx + 0.2])
    ax.set_ylim([-0.2, ny + 0.2])
    ax.set_zlim([-0.2, nz + 0.2])
    ax.set_aspect('equal')
    ax.grid(True)
    ax.legend()
    plt.show()

def plot_fcc_unit_cell():
    """
    FCC의 단일 단위 셀(conventional unit cell)과 원자 공유 개념을 시각화합니다.
    """
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(projection='3d')

    # 1. 꼭짓점 원자 (8개, 각 1/8 기여)
    corners = np.array([
        [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
        [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1]
    ])
    ax.scatter(corners[:, 0], corners[:, 1], corners[:, 2], s=200, c='blue', label='Corners (8 * 1/8 = 1 atom)')

    # 2. 면 중심 원자 (6개, 각 1/2 기여)
    faces = np.array([
        [0.5, 0.5, 0], [0.5, 0.5, 1],
        [0.5, 0, 0.5], [0.5, 1, 0.5],
        [0, 0.5, 0.5], [1, 0.5, 0.5]
    ])
    ax.scatter(faces[:, 0], faces[:, 1], faces[:, 2], s=200, c='red', label='Faces (6 * 1/2 = 3 atoms)')

    # 단위 셀의 경계를 선으로 그립니다.
    for s, e in [[(0,0,0),(1,0,0)], [(0,0,0),(0,1,0)], [(0,0,0),(0,0,1)],
                 [(1,0,0),(1,1,0)], [(1,0,0),(1,0,1)], [(0,1,0),(1,1,0)],
                 [(0,1,0),(0,1,1)], [(0,0,1),(1,0,1)], [(0,0,1),(0,1,1)],
                 [(1,1,0),(1,1,1)], [(1,0,1),(1,1,1)], [(0,1,1),(1,1,1)]]:
        ax.plot3D(*zip(s, e), color="gray", linestyle='--')

    # 축 레이블 및 제목 설정
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('FCC Unit Cell (Total Atoms = 1 + 3 = 4)', fontsize=16)

    # 축의 범위를 설정하여 정육면체 형태로 보이게 합니다.
    ax.set_xlim([-0.2, 1.2])
    ax.set_ylim([-0.2, 1.2])
    ax.set_zlim([-0.2, 1.2])
    ax.set_aspect('equal')
    ax.legend()
    ax.grid(False) # 배경 격자는 끔

    # 플롯을 화면에 표시합니다.
    plt.show()

if __name__ == "__main__":
    # --- 격자 크기 설정 ---
    # 여기에서 x, y, z 방향의 단위 셀 개수를 조절할 수 있습니다.
    NX, NY, NZ = 5, 4, 2

    # 지정된 크기의 FCC 격자 원자 위치를 생성합니다.
    corner_atoms, face_center_atoms = generate_fcc_lattice(nx=NX, ny=NY, nz=NZ)
    
    # 생성된 격자를 시각화합니다.
    plot_title = f'{NX}x{NY}x{NZ} FCC Lattice (Colored by Type)'
    plot_lattice(corners=corner_atoms, face_centers=face_center_atoms, title=plot_title, grid_dims=(NX, NY, NZ))

    # 단일 단위 셀을 보려면 아래 줄의 주석을 해제하고 위 코드를 주석 처리하세요.
    # plot_fcc_unit_cell()

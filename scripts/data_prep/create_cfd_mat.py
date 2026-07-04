from pathlib import Path
import argparse
import re

import numpy as np
from scipy.io import savemat


def get_timestep(path):
    m = re.search(r"t=(\d+)", path.name)
    return int(m.group(1)) if m else -1


def parse_args():
    project_root = Path(__file__).resolve().parents[2]

    parser = argparse.ArgumentParser(
        description="Create a cropped/downsampled MATLAB CFD snapshot file from extracted .txt snapshots."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(r"C:\Users\edric\Desktop\2d_cfd_data"),
        help="Folder containing extracted .txt snapshot files.",
    )
    parser.add_argument(
        "--mat-file",
        type=Path,
        default=project_root / "data" / "processed" / "cylinder_2D_uv_ONE_MIDDLE_Z_cropped_downsampled.mat",
        help="Output .mat file.",
    )
    parser.add_argument("--max-states", type=int, default=100000)
    parser.add_argument("--x-min", type=float, default=4.5)
    parser.add_argument("--x-max", type=float, default=7.5)
    parser.add_argument("--y-min", type=float, default=3.0)
    parser.add_argument("--y-max", type=float, default=5.0)
    parser.add_argument("--downsample-stride", type=int, default=4)
    parser.add_argument("--max-snapshots", type=int, default=None)
    return parser.parse_args()


def main():
    args = parse_args()

    # ======================================================
    # 1. Folder containing extracted .txt snapshot files
    # ======================================================
    data_dir = args.data_dir
    mat_file = args.mat_file
    mat_file.parent.mkdir(parents=True, exist_ok=True)

    # ======================================================
    # 2. User settings
    # ======================================================
    max_states = args.max_states
    max_points = max_states // 2

    x_min, x_max = args.x_min, args.x_max
    y_min, y_max = args.y_min, args.y_max

    downsample_stride = args.downsample_stride

    max_snapshots = args.max_snapshots
    # max_snapshots = 50

    # ======================================================
    # 3. Find and sort files
    # ======================================================
    files = list(data_dir.glob("*.txt"))
    files = sorted(files, key=get_timestep)

    if max_snapshots is not None:
        files = files[:max_snapshots]

    if len(files) == 0:
        raise RuntimeError(f"No .txt files found in {data_dir}")

    print(f"Found {len(files)} snapshot files")

    # ======================================================
    # 4. Read first snapshot and create ONE-Z crop/downsample mask
    # ======================================================
    A0 = np.loadtxt(files[0], comments="#", dtype=np.float32)

    # columns: x y z u v w
    x_all = A0[:, 0]
    y_all = A0[:, 1]
    z_all = A0[:, 2]

    # Pick the middle z-plane automatically
    unique_z = np.unique(z_all)
    z_target = unique_z[len(unique_z) // 2]
    z_tol = 1e-6

    print("\nUsing ONE z-plane only:")
    print(f"z_target = {float(z_target):.6f}")

    crop_mask = (
        (x_all >= x_min)
        & (x_all <= x_max)
        & (y_all >= y_min)
        & (y_all <= y_max)
        & (np.abs(z_all - z_target) <= z_tol)
    )

    idx_crop = np.where(crop_mask)[0]

    print("\nInitial crop:")
    print(f"x: {x_min} to {x_max}")
    print(f"y: {y_min} to {y_max}")
    print(f"z: {float(z_target):.6f}")
    print(f"Cropped points before downsample: {len(idx_crop)}")
    print(f"Cropped states before downsample: {2 * len(idx_crop)}")

    if len(idx_crop) == 0:
        raise RuntimeError("Crop found 0 points. Increase z_tol or check x/y/z ranges.")

    # Downsample inside the single z-plane crop
    idx_keep = idx_crop[::downsample_stride]

    # Optional hard cap in case still too large
    if len(idx_keep) > max_points:
        print(f"\nDownsampled crop has {len(idx_keep)} points.")
        print(f"Trimming to {max_points} points so states <= {max_states}.")
        idx_keep = idx_keep[:max_points]

    x = x_all[idx_keep].astype(np.float32)
    y = y_all[idx_keep].astype(np.float32)
    z = z_all[idx_keep].astype(np.float32)

    npoints = len(idx_keep)
    nt = len(files)

    print("\nFinal selected data:")
    print(f"npoints: {npoints}")
    print(f"state size [u; v]: {2 * npoints}")
    print(f"snapshots: {nt}")
    print(f"estimated X memory GB: {(2 * npoints * nt * 4) / 1e9:.4f}")
    print(f"x range kept: {float(x.min()):.4f} to {float(x.max()):.4f}")
    print(f"y range kept: {float(y.min()):.4f} to {float(y.max()):.4f}")
    print(f"z range kept: {float(z.min()):.4f} to {float(z.max()):.4f}")

    # ======================================================
    # 5. Build X = states x time
    # State vector = [u values; v values]
    # ======================================================
    X = np.empty((2 * npoints, nt), dtype=np.float32)
    tvals = np.empty(nt, dtype=np.int32)

    for k, f in enumerate(files):
        A = np.loadtxt(f, comments="#", dtype=np.float32)

        if A.shape[0] != A0.shape[0]:
            raise RuntimeError(f"Point count changed in {f.name}")

        u = A[idx_keep, 3].astype(np.float32)
        v = A[idx_keep, 4].astype(np.float32)

        X[:, k] = np.concatenate([u, v])
        tvals[k] = get_timestep(f)

        if (k + 1) % 10 == 0 or k == 0:
            print(f"Loaded {k + 1}/{nt}")

    print("\nFinal X shape:", X.shape)
    print("Final X memory GB:", X.nbytes / 1e9)

    # ======================================================
    # 6. Save MATLAB file
    # ======================================================
    savemat(
        mat_file,
        {
            "X": X,
            "x": x,
            "y": y,
            "z": z,
            "tvals": tvals,
            "npoints": npoints,
            "x_min": x_min,
            "x_max": x_max,
            "y_min": y_min,
            "y_max": y_max,
            "z_target": z_target,
            "z_tol": z_tol,
            "downsample_stride": downsample_stride,
            "description": "Cylinder wake ONE middle z-plane velocity snapshots. Cropped and downsampled. X = [u; v], states x time.",
        },
        do_compression=True,
    )

    print("\nSaved MATLAB file to:")
    print(mat_file)


if __name__ == "__main__":
    main()


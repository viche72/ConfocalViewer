#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
lsm_to_npz_v4.py
----------------
Convert Zeiss .LSM (and optionally .TIF/.TIFF) confocal stacks into a browser-friendly
NPZ bundle for the HTML viewer. The NPZ contains:
  - r.npy, g.npy, b.npy  (float32 arrays, shape Z×Y×X, C-order, values in [0, 1])
  - meta.json            (voxel sizes, downsample factor, background sigma, channel map, filename)

Processing pipeline (per channel):
  1) Gaussian background subtraction (σ = --sigma pixels; set 0 to disable)
  2) Robust percentile normalization to [0,1] (0.5–99.5th percentile)
  3) XY downsampling by an integer factor (--ds), Z unchanged

Channel mapping:
  Use --map R,G,B to pick which input channel index becomes R, G, and B.
  Example: --map 2,0,1 means “output Red = source 2; output Green = source 0; output Blue = source 1”.

Typical usage:
  Single file → single NPZ:
    python lsm_to_npz_v4.py "C:\data\LNCaP-2hr.lsm" "C:\out\LNCaP-2hr.npz" --ds 6 --sigma 12 --map 2,0,1

  Batch a folder (LSM only, non-recursive):
    python lsm_to_npz_v4.py "C:\data\stacks" "C:\out\npz" --ds 6 --sigma 12 --map 2,0,1

  Recurse subfolders and include TIFFs as well:
    python lsm_to_npz_v4.py "C:\data\stacks" "C:\out\npz" --recursive --include-tiff --ds 6 --sigma 12

Notes:
  • Voxel sizes are read from LSM metadata (meters → micrometers). If unavailable, falls back to 1.0 µm.
  • The viewer multiplies X/Y voxel sizes by --ds at render time to preserve geometry.
  • This version expects at least 3 channels. If you need 1–2 channel support or time series selection,
    use v5 (which injects missing channels); happy to include that as well.

Author: you :)
License: MIT
"""

import sys, os, json, zipfile, io, argparse
from pathlib import Path
import numpy as np
import tifffile
from scipy.ndimage import gaussian_filter
from numpy.lib.format import write_array


# ---------------------------- Image preprocessing --------------------------------

def normalize(img, lo=0.5, hi=99.5):
    """
    Percentile-based normalization to [0,1].
    lo, hi are percentiles (float). Robust to outliers.
    """
    lo_v, hi_v = np.percentile(img, [lo, hi])
    if hi_v == lo_v:
        return np.zeros_like(img, dtype=np.float32)
    x = (img.astype(np.float32) - lo_v) / (hi_v - lo_v)
    return np.clip(x, 0, 1)


def subtract_bg(img, sigma=12):
    """
    Estimate background via Gaussian blur (sigma in pixels), subtract, and floor at 0.
    Set sigma=0 to disable subtraction.
    """
    if sigma and sigma > 0:
        bg = gaussian_filter(img.astype(np.float32), sigma=sigma)
        out = img.astype(np.float32) - bg
        out[out < 0] = 0
        return out
    return img.astype(np.float32)


# ---------------------------- Core conversion logic ------------------------------

def convert_one(in_path: Path, out_npz: Path, ds=6, sigma=12, ch_map=(2, 0, 1)):
    """
    Convert a single LSM/TIFF to NPZ with r.npy, g.npy, b.npy and meta.json.

    Assumptions:
      - tif.series[0] is the main 4D array with axes Z,C,Y,X (typical for LSM).
      - There are ≥ 3 channels. (For 1–2 channel datasets, consider v5.)
    """

    with tifffile.TiffFile(str(in_path)) as tif:
        series = tif.series[0]
        vol = series.asarray()             # expected shape: (Z, C, Y, X)
        meta = getattr(tif, "lsm_metadata", {}) or {}

    if vol.ndim != 4:
        raise ValueError(f"{in_path.name}: expected 4D (Z,C,Y,X), got {vol.shape}")

    Z, C, Y, X = vol.shape
    if C < 3:
        raise ValueError(f"{in_path.name}: need 3 channels (got {C}).")

    # Physical voxel sizes (meters in LSM → micrometers)
    vx_um = float(meta.get("VoxelSizeX", 1.0)) * 1e6
    vy_um = float(meta.get("VoxelSizeY", 1.0)) * 1e6
    vz_um = float(meta.get("VoxelSizeZ", 1.0)) * 1e6

    def clean(ch):
        """
        Apply background subtraction and normalization per Z slice, then
        downsample in X/Y by factor ds. Return float32 Z×Y×X in [0,1].
        """
        stack = np.stack([subtract_bg(vol[z, ch], sigma) for z in range(Z)], axis=0)
        stack = np.stack([normalize(stack[z]) for z in range(Z)], axis=0)
        return stack[:, ::ds, ::ds].astype('float32', copy=False)

    # Channel mapping: output R,G,B are selected source-channel indices
    src_for_r, src_for_g, src_for_b = ch_map
    r = clean(src_for_r)
    g = clean(src_for_g)
    b = clean(src_for_b)

    # Metadata bundled into NPZ
    meta_out = {
        "vx_um": float(vx_um),
        "vy_um": float(vy_um),
        "vz_um": float(vz_um),
        "ds": int(ds),
        "sigma": float(sigma),
        "file": in_path.name,
        "channel_map": list(map(int, ch_map))
    }

    # Write standard NPZ with NPY v1.0 members (so the browser parser stays simple)
    out_npz.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(str(out_npz), "w", compression=zipfile.ZIP_DEFLATED) as zf:

        def write_npy(arcname, arr):
            bio = io.BytesIO()
            arr_c = np.ascontiguousarray(arr, dtype=np.float32)
            write_array(bio, arr_c, version=(1, 0))  # NPY v1.0 header
            zf.writestr(arcname, bio.getvalue())

        write_npy("r.npy", r)
        write_npy("g.npy", g)
        write_npy("b.npy", b)
        zf.writestr("meta.json", json.dumps(meta_out).encode("utf-8"))

    print(f"[OK] {in_path.name} -> {out_npz.name}  |  Z×Y×X={tuple(r.shape)}  |  "
          f"voxel µm=({vx_um:.3f},{vy_um:.3f},{vz_um:.3f})  |  map={ch_map}")


# ---------------------------- File discovery -------------------------------------

def discover_files(root: Path, include_tiff=False, recursive=False):
    """
    Return a sorted list of input files to process.

    include_tiff:
      False → only .lsm/.LSM
      True  → also include .tif/.tiff/.TIF/.TIFF

    recursive:
      False → only direct children of 'root' if it is a folder
      True  → walk subfolders
    """
    exts = [".lsm", ".LSM"]
    if include_tiff:
        exts += [".tif", ".tiff", ".TIF", ".TIFF"]

    if root.is_file():
        return [root]

    files = []
    if recursive:
        for p in root.rglob("*"):
            if p.is_file() and p.suffix in exts:
                files.append(p)
    else:
        for p in root.iterdir():
            if p.is_file() and p.suffix in exts:
                files.append(p)
    return sorted(files)


# ---------------------------- CLI ------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Convert LSM (and optionally TIFF) into NPZ (r,g,b + meta.json) for the HTML volume viewer."
    )
    ap.add_argument("input",
                    help="Path to a single file OR a folder to process.")
    ap.add_argument("output",
                    help="Output NPZ path (for a single file) OR output folder (for batch).")

    # Ranges and guidance are included in help strings:
    ap.add_argument("--ds", type=int, default=6,
                    help="XY downsample factor (int, recommended 1–16; default 6). "
                         "Higher = faster/smaller, lower = sharper XY at higher memory cost.")

    ap.add_argument("--sigma", type=float, default=12,
                    help="Gaussian background sigma in pixels (float, recommended 0–50; default 12). "
                         "Use 0 to disable background subtraction.")

    ap.add_argument("--map", type=str, default="2,0,1",
                    help="Channel map 'r,g,b' = source indices (three comma-separated non-negative ints). "
                         "Example: 2,0,1 → Red←Ch2, Green←Ch0, Blue←Ch1. "
                         "Each index should be in [0, C-1], where C is # of channels in the file.")

    ap.add_argument("--recursive", action="store_true",
                    help="Recurse into subfolders (only relevant when input is a folder).")

    ap.add_argument("--include-tiff", action="store_true",
                    help="Also include .tif/.tiff files (default processes only .lsm).")

    args = ap.parse_args()

    # Parse and validate --map
    try:
        ch_map = tuple(int(x.strip()) for x in args.map.split(","))
        if len(ch_map) != 3 or any(i < 0 for i in ch_map):
            raise ValueError
    except Exception:
        raise SystemExit("--map must be THREE comma-separated non-negative integers, e.g. '2,0,1'")

    in_path = Path(args.input)
    out_path = Path(args.output)

    files = discover_files(in_path, include_tiff=args.include_tiff, recursive=args.recursive)
    if not files:
        print("No matching files found.")
        sys.exit(1)

    if in_path.is_file():
        # Single input file: output may be a file (.npz) or a folder
        if out_path.is_dir() or out_path.suffix.lower() != ".npz":
            out_npz = out_path / (in_path.stem + ".npz")
        else:
            out_npz = out_path
        convert_one(in_path, out_npz, ds=args.ds, sigma=args.sigma, ch_map=ch_map)
    else:
        # Batch mode
        out_path.mkdir(parents=True, exist_ok=True)
        for p in files:
            out_npz = out_path / (p.stem + ".npz")
            try:
                convert_one(p, out_npz, ds=args.ds, sigma=args.sigma, ch_map=ch_map)
            except Exception as e:
                print(f"[ERR] {p.name}: {e}")


if __name__ == "__main__":
    main()

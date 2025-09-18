
---

````markdown
# Zeiss LSM → Browser-based 3-D Volume Viewer

This repository provides a **complete workflow** to transform Zeiss confocal
`.lsm` stacks (and optionally standard `.tif/.tiff` stacks) into an
**interactive 3-D volume** that you can explore directly in your browser—no server required.

---

## Contents

| File | Purpose |
|------|--------|
| **`lsm_to_npz_v4.py`** | Python converter: processes Zeiss LSM or TIFF stacks and produces a compact **`.npz` bundle** (`r.npy`, `g.npy`, `b.npy`, `meta.json`) for the viewer. |
| **`volume_viewer.html`** | Stand-alone HTML/JavaScript viewer: load the `.npz` files locally and interactively explore the 3-D data in a WebGL browser window with per-channel color, contrast/gamma controls, orthogonal MIPs, and PNG snapshot export. |

---

## Quick Start

### 1. Install prerequisites
Python ≥3.9 recommended.

```bash
pip install numpy scipy tifffile
````

*No extra dependencies are needed for the HTML viewer; it runs entirely in the browser.*

---

### 2. Convert LSM/TIFF to NPZ

```bash
python lsm_to_npz.py INPUT OUTPUT [options]
```

**Examples**

```bash
# Single file → single NPZ (explicit file path)
python lsm_to_npz.py "C:\data\LNCaP-2hr.lsm" "C:\out\LNCaP-2hr.npz" --ds 6 --sigma 12 --map 2,0,1

# Folder (LSM only, non-recursive) → NPZs with same base name
python lsm_to_npz.py "C:\data\stacks" "C:\out\npz" --ds 6 --sigma 12 --map 2,0,1

# Recurse subfolders and include standard TIFFs
python lsm_to_npz.py "C:\data\stacks" "C:\out\npz" --recursive --include-tiff --ds 6 --sigma 12
```

---

### 3. View interactively

1. **Double-click** `volume_viewer.html` (Chrome/Edge with hardware acceleration enabled).
2. Click **Load NPZ(s)** and select one or more `.npz` files.
3. Use the **Dataset** dropdown to switch between loaded stacks.

---

## `lsm_to_npz.py`

### What it does

* Reads Zeiss `.lsm` (and optionally `.tif/.tiff`) via **tifffile**.
* Extracts voxel sizes from LSM metadata (meters → **µm**).
* For each output channel (**R,G,B**):

  1. **Background subtraction** with a Gaussian blur (`--sigma` px).
  2. **Percentile normalization** to \[0,1] (0.5–99.5th percentile).
  3. **XY downsampling** by integer factor (`--ds`), Z unchanged.
* Writes a single `.npz` containing:

  * `r.npy`, `g.npy`, `b.npy` — float32 arrays of shape **Z × Y × X**.
  * `meta.json` — voxel sizes (`vx_um`, `vy_um`, `vz_um`), `ds`, `sigma`, `channel_map`, and input filename.

> **Requirement:** The first series of the LSM must have at least **3 channels**.
> (If you need support for 1–2 channel datasets or time-lapse data, use the v5 script described in the documentation.)

---

### Command-line options

| Flag             | Type         | Description                                                                                 | Recommended Range  | Default |
| ---------------- | ------------ | ------------------------------------------------------------------------------------------- | ------------------ | ------- |
| `input`          | path         | Single LSM/TIFF file or folder to process.                                                  | —                  | —       |
| `output`         | path         | Output `.npz` file (if input is a single file) **or** output folder (if input is a folder). | —                  | —       |
| `--ds`           | int          | **XY downsample factor**. Larger → smaller/faster files, smaller → higher XY resolution.    | **1–16**           | 6       |
| `--sigma`        | float        | **Gaussian σ** for background subtraction (pixels). 0 disables subtraction.                 | **0–50**           | 12      |
| `--map`          | r,g,b (ints) | **Channel mapping**: which input channels become Red, Green, Blue. Indices in `[0, C-1]`.   | depends on dataset | `2,0,1` |
| `--recursive`    | flag         | Recurse into subfolders when input is a folder.                                             | —                  | off     |
| `--include-tiff` | flag         | Also process `.tif/.tiff` files.                                                            | —                  | off     |

---

### Tips for increased fidelity

* **Lower `--ds`** to 4, 2 or even 1 to capture more XY detail (file size and GPU memory use grow \~quadratically).
* Keep `--sigma` modest (8–12 px) to avoid blurring fine structures.
* Crop regions of interest before conversion if you only need a subset of the field.

---

## `volume_viewer_client_v3_2.html`

Open this HTML file directly in a modern browser—no server needed.

### Features

* **Multi-dataset loading** with dataset dropdown.
* **Per-channel controls**

  * Visibility toggle
  * Color picker
  * `isomin` (threshold) and `opacity` sliders
* **Global tone controls**

  * `Lo`, `Hi`, `Gamma` (applied to all channels):

    ```
    v' = clip((v − Lo)/(Hi − Lo), 0..1)^Gamma
    ```
* **Orthogonal maximum-intensity projections** (XY, XZ, YZ).
* **PNG snapshot** export (1200×800).

### Recommended UI ranges

| Control     | Purpose                           | Typical Range        |
| ----------- | --------------------------------- | -------------------- |
| Isosurfaces | number of isosurfaces per channel | 5–30                 |
| Isomin      | isosurface threshold              | 0.05–0.5             |
| Opacity     | per-channel alpha                 | 0.02–0.15            |
| Lo / Hi     | global contrast window            | Lo 0–0.5, Hi 0.5–1.0 |
| Gamma       | mid-tone boost                    | 0.5–3.0              |

### Geometry

* Uses **true physical voxel sizes** from `meta.json`.
* Multiplies X/Y voxel sizes by `ds` (the downsample factor used during conversion) so the 3-D aspect ratio remains correct.

---

## Performance & Fidelity

* **More detail** → lower `--ds`, higher `surface_count`, tighter Lo/Hi window, Gamma \~1.3–1.6.
* **Faster interaction** → higher `--ds`, lower `surface_count`, raise `Isomin`, toggle off unused channels.
* A discrete GPU and Chrome/Edge with hardware acceleration are recommended when `--ds ≤ 2` on 2k×2k stacks.

---

## Troubleshooting

* **File not found**: Use your actual local file paths (e.g. `C:\Users\...`), not `/mnt/data/...` which is only for cloud notebooks.
* **Volume looks stretched in Z**: Check the *meta pill* in the viewer; X/Y voxel sizes should be `vx_um * ds`, Z = `vz_um`.
* **No color**: Make sure `r.npy`, `g.npy`, `b.npy` exist and adjust per-channel `isomin`/`opacity`.

---

## License

MIT – use, modify, and distribute freely. Please credit the underlying libraries:

* [tifffile](https://pypi.org/project/tifffile/) for LSM reading
* [SciPy](https://scipy.org/) for Gaussian filtering
* [Plotly](https://plotly.com/javascript/) and [fflate](https://github.com/101arrowz/fflate) for the browser viewer.

---

## Citation

If these tools contribute to a publication or presentation, please cite this repository and the libraries listed above.

---

## Roadmap / Advanced Options

* **High-fidelity rendering**: lower `--ds` to 1 or 2 and increase `surface_count` to 20–30.
* **Large datasets**: consider cropping or using a multiscale format such as [OME-Zarr](https://ome-zarr.readthedocs.io/) with [vizarr](https://vizarr.org/).
* **Mesh export**: for publication-quality figures, add a Python step with `skimage.measure.marching_cubes` to generate `.ply` or `.glb` isosurface meshes.

---

*Happy exploring your confocal stacks in 3-D!*



# Workflow

This document summarizes the recommended workflow for lab members reusing this repository.

## 1. Prepare the Environment

Use the unified Python 3.13.13 environment:

```bash
conda env create -f tests/environment.yaml
```

This merges both former environment files and pins Suite2p 1.1.0. OpenCV uses
`opencv-python-headless` to share the same package as Cellpose; project interaction
and the Suite2p GUI use PyQt6. Do not additionally install `opencv-python`, because
both distributions write the same `cv2` files.

Validation commands (from the repository root after activation):

```bash
python -m pip check
python tests/validate_environment.py
python -m unittest discover -s tests -v
python tests/smoke_suite2p.py
```

The smoke test uses temporary synthetic images, runs real CPU registration, and
saves its report under `reports/environment_validation/`. Installed versions are
recorded in `reports/environment_validation/requirements-lock.txt`; use this as a
pip constraints file with the environment's declared requirements when reproducing
the tested package set. NWB import/export is disabled and its optional `pynwb`
dependency is not included.

Use `ID-mRSFC-ISOI_suite2p` for preprocessing:

```bash
conda activate ID-mRSFC-ISOI_suite2p
python scripts/01_preprocess.py -c configs/default.yaml
```

Use the same environment for downstream analysis:

```bash
conda activate ID-mRSFC-ISOI_suite2p
python scripts/02_correlation.py -c configs/default.yaml
python scripts/03_identification.py -c configs/default.yaml
python scripts/04_statistics.py -c configs/default.yaml
```

## 2. Prepare Raw Data

Raw data should follow the config-defined layout:

```text
<raw>/<dataset>/<trial>/<session>/<mouse>/
|-- info.txt
|-- img*.bin
`-- ai*.bin
```

With `configs/default.yaml`, this means:

```text
data/raw/WT/Trial1/Session1/A01/
```

With `configs/OIS.yaml` (`dataset: ""`, trials `00w` through `10w`), this means:

```text
F:/OIS_Awake_Mice/00w/Run01/A01/
```

Only the raw input layout uses `trial/session/mouse`. Generated per-recording outputs use `session/mouse/trial`, as expected by the correlation stage.

## 3. Prepare Masks

ROI masks are loaded from `paths.mask`, usually `data/mask/`.

Each mask should be a PNG image where nonzero pixels define one unilateral ROI. The code mirrors each mask left/right when `extract_roi(..., mode="symmetric")` is used.

The shared chamber mask is stored alongside the ROI masks:

```text
data/mask/Chamber_Mask.png
```

Filename case differences (e.g. `Chamber_mask.png`) are accepted. It is used for all recordings and excluded from ROI extraction. Its dimensions must match the aligned frames. All datasets using the same `paths.mask` share this mask.

## 4. Run Preprocessing

Suite2p writes directly under each processed recording folder as `suite2p/plane0/reg_tif/`. It receives an explicit input file list containing only `<mouse>_<trial>.tif`, so re-running does not include `_merged.tif`. Existing outputs in the old `OIS_Corrected` directory are not migrated or deleted.

Preprocessing performs:

1. Raw OIS binary to TIFF conversion.
2. Suite2p motion correction.
3. Registered TIFF merging.
4. Interactive atlas alignment.
5. dR/R, detrending, global signal regression, normalization, filtering, and downsampling.

Expected output per trial:

```text
data/processed/<dataset>/<session>/<mouse>/<trial>/
|-- suite2p/                 # registration outputs for this session/mouse/trial
|-- <mouse>_<trial>.tif
|-- <mouse>_<trial>_merged.tif
|-- <mouse>_<trial>_Affine.npy
|-- <mouse>_<trial>_Average_After.png
`-- <mouse>_<trial>_FiltDRR.npy
```

Atlas alignment is interactive:

1. Select Bregma.
2. Select Lambda.
3. Select P3.
4. Select P4.
5. Select P5.
6. Press Enter in the figure window.
7. Type `y` in the terminal to accept, or `n` to reselect.

## 5. Run FC Matrix Generation

`02_correlation.py` reads:

```text
data/processed/<dataset>/<session>/<mouse>/<trial>/<mouse>_<trial>_FiltDRR.npy
```

and writes:

```text
data/processed/<dataset>/FC_Matrix/<session>_<mouse>_<trial>_FC_Matrix_<atlas_tag>_excluded.npy
data/processed/<dataset>/<session>_correlation_matrices_<atlas_tag>_excluded.pkl
data/output/<dataset>/FC_Matrix/<session>_<mouse>_<trial>_FC_Matrix_<atlas_tag>_excluded.png
```

## 6. Run Identification and Statistics

`03_identification.py` loads session-level pkl files according to `Session_ids` in the active config. The config order is preserved.

The script writes similarity, edge-contribution, and summary statistic outputs under:

```text
data/output/<dataset>/
```

Then run `04_statistics.py` to summarize `statistics.pkl`.

## 7. Run Supplementary Structure Analysis

`05_supplement.py` is optional and should be run after `03_identification.py` has produced:

```text
data/output/<dataset>/statistics.pkl
```

Run it in the main analysis environment:

```bash
conda activate ID-mRSFC-ISOI_suite2p
python scripts/05_supplement.py -c configs/default.yaml
```

The script builds a structure-space summary from the identification results:

- DD: diagonal dominance of the subject similarity matrix.
- ND: normalized distance to the diagonal.
- RM-ANOVA: compares DD and ND across configured groups.
- MANOVA: jointly compares DD and ND across configured groups.

Before running, confirm these config fields:

- `supplement.statistics_file`: statistics file name under `data/output/<dataset>/`.
- `supplement.structure_atlases`: atlas conditions to include, using the same tags generated from `variables.Atlas`.
- `supplement.structure_group_names`: display names paired one-to-one with `structure_atlases`.
- `supplement.exclude`: rows to ignore, such as `AUDv`.
- `visualization.supplement_dd_box`, `visualization.supplement_nd_box`, and `visualization.supplement_structure_space`: plot styles.

Expected outputs are written under:

```text
data/output/<dataset>/Stat/
|-- boxplot_DD.png
|-- boxplot_ND.png
`-- structure_space.png
```

Statistical tables are written to the active log file.

## 8. Practical Notes

- Run commands from the repository root.
- Keep the exact config used for a run if you need reproducibility.
- Preprocessing can require substantial disk space because Suite2p writes binary and registered TIFF outputs.
- Do not commit generated files under `data/processed`, `data/output`, `data/raw`, or `logs`.


## Cumulative analysis and preprocessing resume

The current pipeline defaults to random session grouping with a reproducible seed.
All analysis output now lives under a seed/mode/configuration-specific directory inside `paths.output`; paths shown earlier in this document describe filenames and stage contracts, with this additional directory inserted before the dataset.

03 saves separate `1Trial/`, `2Trial/`, ... results. 1Trial has valid identification results but no within-session split-half Edge estimate; its Edge result records this explicitly. Edge estimates at N>=2 use only the first N configured trials.

01 resumes at the first incomplete or changed step, validates its new outputs, then reruns dependent downstream steps. A complete step 4 is skipped without opening the interactive alignment UI. Details and limitations are in [config.md](config.md#analysis-output-and-resume).

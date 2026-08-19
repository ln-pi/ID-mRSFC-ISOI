# Workflow

This document summarizes the recommended workflow for lab members reusing this repository.

## 1. Prepare Environments

Use two conda environments:

```bash
conda env create -f neuro_suite2p_env.yaml
conda env create -f environment.yaml
```

Use `neuro_suite2p_env` only for preprocessing:

```bash
conda activate neuro_suite2p_env
python scripts/01_preprocess.py -c configs/default.yaml
```

Use `ID-mRSFC-ISOI` for downstream analysis:

```bash
conda activate ID-mRSFC-ISOI
python scripts/02_correlation.py -c configs/default.yaml
python scripts/03_identification.py -c configs/default.yaml
python scripts/04_statistics.py -c configs/default.yaml
```

## 2. Prepare Raw Data

Raw data should follow the config-defined layout:

```text
<raw>/<dataset>/<session>/<mouse>/<trial>/
|-- info.txt
|-- img*.bin
`-- ai*.bin
```

With `configs/default.yaml`, this means:

```text
data/raw/WT/Session1/A01/Trial1/
```

With `configs/OIS.yaml`, this means:

```text
F:/OIS_Awake_Mice/WT/Run01/41#/Trial1/
```

## 3. Prepare Masks

ROI masks are loaded from `paths.mask`, usually `data/mask/`.

Each mask should be a PNG image where nonzero pixels define one unilateral ROI. The code mirrors each mask left/right when `extract_roi(..., mode="symmetric")` is used.

Each processed trial folder must also contain:

```text
Chamber_Mask.png
```

before `drr_narrow()` can run.

## 4. Run Preprocessing

Preprocessing performs:

1. Raw OIS binary to TIFF conversion.
2. Suite2p motion correction.
3. Registered TIFF merging.
4. Interactive atlas alignment.
5. dR/R, detrending, global signal regression, normalization, filtering, and downsampling.

Expected output per trial:

```text
data/processed/<dataset>/<session>/<mouse>/<trial>/
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
conda activate ID-mRSFC-ISOI
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

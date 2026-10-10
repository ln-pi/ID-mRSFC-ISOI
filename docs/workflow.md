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
<raw>/<dataset>/<week>/<run>/<mouse>/
|-- info.txt
|-- img*.bin
`-- ai*.bin
```

With `configs/default.yaml`, this means:

```text
data/raw/WT/00w/Run01/A01/
```

With `configs/OIS.yaml` (`dataset: ""`, configured weeks and runs), this means:

```text
F:/OIS_Awake_Mice/00w/Run01/N01/
```

Raw and processed inputs may use any permutation of week/mouse/run directory levels. The configured IDs resolve each recording. Duplicate sources are reported as conflicts. New preprocessing outputs use week/mouse/run; a single existing recording folder can be resumed.

## 3. Prepare Masks

ROI masks are loaded from `paths.mask`, usually `data/mask/`.

Each mask should be a PNG image where nonzero pixels define one unilateral ROI. The code mirrors each mask left/right when `extract_roi(..., mode="symmetric")` is used.

The shared chamber mask is stored alongside the ROI masks:

```text
data/mask/Chamber_Mask.png
```

Filename case differences (e.g. `Chamber_mask.png`) are accepted. It is used for all recordings and excluded from ROI extraction. Its dimensions must match the aligned frames. All datasets using the same `paths.mask` share this mask.

## 4. Run Preprocessing

Suite2p writes directly under each processed recording folder as `suite2p/plane0/reg_tif/`. It receives an explicit input file list containing only `<mouse>_<run>.tif`, so re-running does not include `_merged.tif`. Existing outputs in the old `OIS_Corrected` directory are not migrated or deleted.

Preprocessing performs:

1. Raw OIS binary to TIFF conversion.
2. Suite2p motion correction.
3. Registered TIFF merging.
4. Interactive atlas alignment.
5. dR/R, detrending, global signal regression, normalization, filtering, and downsampling.

Expected output per recording:

```text
data/processed/<dataset>/<week>/<mouse>/<run>/
|-- suite2p/                 # registration outputs for this week/mouse/run
|-- <mouse>_<run>.tif
|-- <mouse>_<run>_merged.tif
|-- <mouse>_<run>_Affine.npy
|-- <mouse>_<run>_Average_After.png
`-- <mouse>_<run>_<week>_FiltDRR.npy
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
data/processed/<dataset>/<week>/<mouse>/<run>/<mouse>_<run>_<week>_FiltDRR.npy
```

and writes:

```text
data/output/<analysis_tag>/<dataset>/Correlation/FC_Matrix/<week>_<mouse>_<run>_FC_Matrix_<atlas_tag>_excluded.npy
data/output/<analysis_tag>/<dataset>/Correlation/<week>_correlation_matrices_<atlas_tag>_excluded.pkl
data/output/<analysis_tag>/<dataset>/Correlation/Figures/<week>_<mouse>_<run>_FC_Matrix_<atlas_tag>_excluded.png
```

## 6. Run Identification and Statistics

`03_identification.py` loads week-level PKLs (mouse rows, run columns).
Set `analysis.split_axis: run|week` and `analysis.selection_mode: single|cumulative`.
Run splits compare Run groups across individual/cumulative weeks; week splits
compare week groups across individual/cumulative runs. Strict/cover always selects
one fixed mouse cohort and week interval before analysis. See docs/config.md.
Edge uses the same split axis and actual groups as ID.

The script writes similarity, edge-contribution, and summary statistic outputs under:

```text
data/output/<analysis_tag>/<dataset>/Identification/
```

Then run `04_statistics.py` to summarize `statistics.pkl`.

## 7. Per-mouse identifiability and optional legacy analysis

The formal pipeline ends at `04_statistics.py`. It reads the full similarity
matrices saved by 03, verifies cohort/analysis metadata and consistency with
`statistics.pkl`, and exports:

- `Statistics/Tables/individual_identifiability.csv`: one row per mouse, atlas and analysis condition.
- `Statistics/individual_identifiability_metadata.json`: definitions, scale, original condition selections and limitations.
- `Statistics/IndividualIdentifiability/`: per-mouse plots of self-match advantage and percentile.

Self-match advantage is the diagonal similarity minus the mean of that mouse's
row/column nonself similarities. Self-match percentile is the average fraction
of other candidates beaten in both directions; exact ties count as half a win.
It is stored on a 0–1 scale (0.5 = 50%). Negative advantage is retained.
Neither quantity is a disease severity score or a new inferential test. Shared
references and candidate composition still matter.

The complete `Identification/<condition>/Numeric/identification_<atlas>.pkl`
files are required: the aggregate alone cannot recover matching percentiles.
Missing files or mismatched cohorts fail explicitly instead of fabricating ranks.

DD/ND has moved out of scripts into `srcs/utils/legacy_structure.py`:

```python
from srcs.utils.config import load_context
from srcs.utils.legacy_structure import run_legacy_structure_analysis

ctx = load_context("configs/default.yaml")
result = run_legacy_structure_analysis(ctx)
```

This explicit function uses existing `supplement` configuration and writes to
`LegacyStructure/Tables/` and `LegacyStructure/Figures/`. Importing it performs
no analysis. Legacy condition-level inference is disabled by default; the
explicit `run_inference=True` option preserves the old calculations solely for
historical comparison, not animal-level disease inference. It is never called
by the formal 01–04 pipeline. Existing historical output folders are not deleted.

## 8. Practical Notes

- Run commands from the repository root.
- Keep the exact config used for a run if you need reproducibility.
- Preprocessing can require substantial disk space because Suite2p writes binary and registered TIFF outputs.
- Do not commit generated files under `data/processed`, `data/output`, `data/raw`, or `logs`.


## Cumulative analysis and preprocessing resume

The pipeline uses a reproducible seed to group analysis.split_axis entries.
All analysis output now lives under a seed/mode/configuration-specific directory inside `paths.output`; paths shown earlier in this document describe filenames and stage contracts, with this additional directory inserted before the dataset.

03 saves separate week/run single conditions or cumulative `NWeek/` or `NRun/` conditions. Edge uses the same axis groups as ID; a single opposite-axis sample can still have a valid Edge estimate.

01 resumes at the first incomplete or changed step, validates its new outputs, then reruns dependent downstream steps. A complete step 4 is skipped without opening the interactive alignment UI. Details and limitations are in [config.md](config.md#preprocessing-resume).

## Missing inputs in stages 01 and 02

Stages 01 and 02 log missing inputs and per-record failures, then continue with
later records. Each run writes a unique `BatchReports/<stage>_<timestamp>_<id>.json`
under the configured analysis output directory, updating it after each record.
The report includes success/missing/failed counts, week, mouse, run, input
path, and error reason (also atlas for stage 02). Stage 01 success includes records
whose complete outputs were reused by the existing resume mechanism.

Stage 02 saves each week/atlas table with `None` in missing or failed cells;
its `complete` attribute indicates whether every configured cell succeeded.
Non-finite FC matrices count as failures. Existing individual output files from
older runs are not deleted; the current table and batch report determine validity.
03 selects a fixed complete mouse cohort using strict/cover and verifies both
the latest complete-scope report and PKL cells. 04/05 use its selected aggregate. Global setup or summary-write errors
still stop the run; user interruption is not swallowed.

## Stage output directories

All stages share the same configuration-specific root:

```text
output/<analysis_tag>/<dataset>/
  Correlation/
    FC_Matrix/                 # 02 numerical arrays and validated caches
    Figures/                   # 02 heatmaps and plot checkpoints
    <week>_correlation_matrices_<atlas>_excluded.pkl
  Identification/
    statistics.pkl             # 03 aggregate; filename configured by analysis.statistics_file
    <condition>/               # week/run ID, NWeek or NRun
    ...
  Statistics/                  # shared 04 and 05 figures
    Tables/                    # 04 statistical tables
      Supplement/              # 05 statistical tables
  BatchReports/
  analysis_config.json
```

04 and 05 read the aggregate from Identification, then write to Statistics.
Existing files in the former output locations are not automatically moved or read.
02 numerical caches remain in Correlation/FC_Matrix; missing Figures are regenerated
from valid caches without repeating numerical computation.

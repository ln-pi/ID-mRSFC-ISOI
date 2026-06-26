# ID-mRSFC-ISOI

This repository contains a Python workflow for mouse resting-state functional connectivity analysis from intrinsic signal optical imaging (ISOI/OIS) data. The pipeline converts raw OIS binary recordings to TIFF stacks, performs motion correction and atlas registration, extracts ROI-level time series, computes functional connectivity (FC) matrices, and evaluates individual-level network identifiability from mRSFC patterns.

The code is organized as research scripts rather than a packaged command-line application. Most parameters are controlled by YAML configuration files in `configs/`.

For lab-member setup and execution, see `docs/workflow.md`. For configuration details, see `docs/config.md`.

## Main Features

- Read OIS200-style raw binary files (`img*.bin`, `ai*.bin`, `info.txt`).
- Convert raw image streams to TIFF stacks.
- Run Suite2p registration in an isolated Python process.
- Interactively align imaging data to the built-in brain atlas.
- Compute dR/R, detrending, global signal regression, band-pass filtering, and downsampling.
- Extract atlas/ROI mean time series from mask images.
- Compute ROI-wise FC matrices.
- Run individual identification analysis from split session/trial FC matrices.
- Estimate edgewise contribution to group consistency and differential power.
- Generate heatmaps, paired plots, boxplots, circular connectivity graphs, and statistical summaries.

## Repository Structure

```text
ID-mRSFC-ISOI/
|-- configs/
|   |-- default.yaml          # Example configuration with Session1-Session4
|   `-- OIS.yaml              # Default config used by scripts unless -c is passed
|-- data/
|   |-- mask                  # Masks for atlas
|   |-- raw                   # Raw OIS data
|   |-- processed             # Processed data: FC matrices and others
|   `-- output                # Outputs: Figures and others
|-- logs/
|-- scripts/
|   |-- 00_main.py            # Placeholder; currently empty
|   |-- 01_preprocess.py      # Raw bin -> TIFF -> Suite2p -> atlas alignment -> FiltDRR
|   |-- 02_correlation.py     # ROI extraction and FC matrix generation
|   |-- 03_identification.py  # Individual identification and edge contribution
|   |-- 04_statistics.py      # Repeated-measures ANOVA and post-hoc tests
|   `-- 05_supplement.py      # Supplementary structure-space analysis
|-- srcs/
|   |-- atlas/                # Built-in atlas assets and masks
|   |-- comp/                 # Alignment, filtering, FC, ID, and edge-contribution tools
|   |-- io/                   # Raw binary readers and data writers
|   |-- prep/                 # Preprocessing and Suite2p wrapper
|   |-- utils/                # Config, ROI extraction, matrix utilities, logger
|   `-- vis/                  # Plotting and circular graph utilities
|-- notebooks/
|   `-- exp.ipynb             # Exploratory notebook
|-- docs/
|   |-- workflow.md           # Lab-member execution workflow
|   `-- config.md             # Configuration field reference
|-- environment.yaml          # Conda environment definition
|-- neuro_suite2p_env.yaml    # Preprocessing environment with suite2p==1.0.0.1
|-- CITATION.cff              # Repository citation metadata
|-- LICENSE                   # MIT license
`-- README.md
```

Generated folders such as `data/` and `logs/` are ignored by Git.

## Environment Setup

This project uses two environments:

- `neuro_suite2p_env`: used for `scripts/01_preprocess.py`, pinned to `suite2p==1.0.0.1` and Python 3.11.
- `ID-mRSFC-ISOI`: used for downstream FC, identification, statistics, plotting, and notebooks.

Create the preprocessing environment:

```bash
conda env create -f neuro_suite2p_env.yaml
conda activate neuro_suite2p_env
```

Create the main analysis environment:

```bash
conda env create -f environment.yaml
conda activate ID-mRSFC-ISOI
```

The preprocessing environment contains the packages needed by `01_preprocess.py` and its direct call chain:

- `suite2p==1.0.0.1` for motion correction.
- `python=3.11`, following Suite2p's current upstream installation recommendation.
- `numpy>=1.24.3`, `scipy`, `pandas`, `tifffile`, and `opencv` for binary/TIFF processing, filtering, masks, and image I/O.
- `matplotlib` and `pyqt` for interactive atlas registration.
- `pyyaml` and `tqdm` for configuration and progress reporting.

The NumPy requirement follows the `suite2p==1.0.0.1` package metadata, which specifies a lower bound and no upper bound.

The main analysis environment is defined by `environment.yaml`. It includes most scientific, graph/network, imaging, notebook, and statistics packages used by `02_correlation.py` through `05_supplement.py`. If your environment does not already contain them, install the shared runtime helpers as well:

```bash
pip install pyyaml tqdm seaborn
```

## Data Layout

The preprocessing script expects raw data to be organized as:

```text
<raw>/<dataset>/<session>/<mouse>/<trial>/
|-- info.txt
|-- img*.bin
`-- ai*.bin
```

For example, using `configs/default.yaml`:

```text
data\raw\WT\Session1\A01\Trial1\
|-- info.txt
|-- img0001.bin
|-- img0002.bin
`-- ...
```

The ROI extraction step expects atlas mask images in `paths.mask`, usually:

```text
data/mask/
|-- <ROI_NAME>.png
|-- <ROI_NAME>.png
`-- ...
```

Each mask is read as a unilateral mask. The code automatically creates left/right ROI pairs by mirroring the mask image.

## Configuration

Scripts read `configs/OIS.yaml` by default. You can pass another config with `-c` or `--config`; the examples below use `configs/default.yaml`.

Important fields, shown with values from `configs/default.yaml`:

```yaml
paths:
  raw: "data\\raw\\"
  processed: "data\\processed\\"
  output: "data\\output\\"
  mask: "data\\mask\\"
  log: "logs\\default.log"

variables:
  dataset: "WT"
  Session_ids: ["Session1", "Session2", "Session3", "Session4"]
  Mouse_ids: ["A01", "A02", "A03", "A04", "A05", "A06", "A07", "A08", "A09", "A10"]
  Trial_ids: ["Trial1", "Trial2", "Trial3", "Trial4", "Trial5", "Trial6", "Trial7"]
  Atlas: [["SM", "VIS", "AUD"], ["PF", "LAT", "MED"], ["None"], ["AUDv"]]
```

`Atlas` entries are exclusion sets used during ROI extraction and identification analysis. For example, `["SM", "VIS", "AUD"]` excludes ROIs whose names contain those labels.

Analysis parameters:

- `ID.group_mode`: defines how four sessions are split into two halves in `build_half_dataframe`.
  - `12`: sessions 1+2 vs 3+4
  - `13`: sessions 1+3 vs 2+4
  - `14`: sessions 1+4 vs 2+3
- `ID.shuffle`: enables permutation testing.
- `ID.n_iter`: shuffle/permutation iteration count used by `03_identification.py`.
- `Edge.repeat`: number of repeated random trial splits for edge contribution.
- `Edge.threshold`: percentile threshold used in circular graph visualization.
- `Edge.cross_only`: whether circular graphs keep cross-module edges only.

Plot appearance is controlled under `visualization` in the active YAML config. This includes FC heatmaps, ID shuffle histograms, paired similarity plots, similarity matrices, edge-contribution circular graphs, boxplots, and statistics summary curves.

## Running the Pipeline

Run commands from the repository root.

### 1. Preprocess Raw Data

Activate the Suite2p environment first:

```bash
conda activate neuro_suite2p_env
```

```bash
python scripts/01_preprocess.py -c configs/default.yaml
```

For each session, mouse, and trial, this script:

1. Reads `img*.bin` and `info.txt`.
2. Saves a raw TIFF stack.
3. Runs Suite2p motion correction.
4. Merges Suite2p registered TIFF chunks.
5. Opens an interactive atlas-alignment window.
6. Saves affine-aligned stacks.
7. Computes filtered dR/R data.

The preprocessing progress bar advances once per completed session/mouse/trial item.

The relevant function chain is:

- `scripts/01_preprocess.py`: loops through `Session_ids`, `Mouse_ids`, and `Trial_ids` from the config.
- `Bin2Tiff()`: reads `img*.bin` with metadata from `info.txt`, merges channel frames, sanitizes invalid frames, and writes `<mouse>_<trial>.tif`.
- `run_isolated_process()`: launches Suite2p in a subprocess using the same active Python interpreter, runs registration with detection/deconvolution disabled, and exports registered TIFF chunks when Suite2p only writes binary output.
- `merge_tiffs()`: collects Suite2p `plane0/reg_tif/*.tif` chunks and writes `<mouse>_<trial>_merged.tif`.
- `Affine_Atlas()`: interactively aligns the average registered frame to the built-in atlas and saves `<mouse>_<trial>_Affine.npy`.
- `drr_narrow()`: applies `Chamber_Mask.png`, computes dR/R, detrends, regresses global signal, normalizes, band-pass filters, downsamples from 30 Hz to 5 Hz, and saves `<mouse>_<trial>_FiltDRR.npy`.

Important: atlas alignment is interactive. You must select five anchors in order:

1. Bregma
2. Lambda
3. P3
4. P4
5. P5

Then press Enter to confirm. After preview, type `y` to accept or `n` to reselect.

Expected outputs include:

```text
data/processed/<dataset>/<session>/<mouse>/<trial>/
|-- <mouse>_<trial>.tif
|-- <mouse>_<trial>_merged.tif
|-- <mouse>_<trial>_Affine.npy
|-- <mouse>_<trial>_Average_After.png
`-- <mouse>_<trial>_FiltDRR.npy
```

### 2. Compute FC Matrices

Switch to the main analysis environment:

```bash
conda activate ID-mRSFC-ISOI
```

```bash
python scripts/02_correlation.py -c configs/default.yaml
```

This script:

1. Loads each filtered dR/R stack.
2. Extracts ROI mean time series using masks from `paths.mask`.
3. Computes Pearson correlation matrices.
4. Saves FC matrices and heatmap figures.

Expected outputs:

```text
data/processed/<dataset>/FC_Matrix/
`-- <session>_<mouse>_<trial>_FC_Maticx_<atlas_tag>_excluded.npy

data/processed/<dataset>/
`-- <session>_correlation_matrices_<atlas_tag>_excluded.pkl

data/output/<dataset>/FC_Matrix/
`-- <session>_<mouse>_<trial>_FC_Matrix_<atlas_tag>_excluded.png
```

### 3. Run Individual Identification

```bash
python scripts/03_identification.py -c configs/default.yaml
```

This script loads packaged FC matrices, splits data into two halves, and computes:

- subject-by-subject similarity matrices,
- individual similarity (`Indiv`),
- group similarity (`Group`),
- relative effect,
- Cohen's d,
- identification accuracy,
- shuffle/permutation results,
- edgewise group consistency (`Phi`),
- edgewise differential power (`DP`).

Correlation matrix package files are loaded from `Session_ids` in the active config, preserving config order.

Expected outputs:

```text
data/output/<dataset>/Similarity/
data/output/<dataset>/EdgeContri/
data/output/<dataset>/Stat/
data/output/<dataset>/statistics.pkl
```

### 4. Run Statistics

```bash
python scripts/04_statistics.py -c configs/default.yaml
```

This script loads:

```text
data/output/<dataset>/statistics.pkl
```

It then generates group/individual trend plots and logs repeated-measures ANOVA plus post-hoc paired tests with FDR correction.

### 5. Supplementary Analysis

```bash
python scripts/05_supplement.py
```

This script defines `StructureAnalyzer`, which compares diagonal dominance (DD) and normalized distance to diagonal (ND) across systems using RM-ANOVA, MANOVA, and structure-space plots.

It reads `statistics.pkl`, compared atlas groups, excluded rows, and plot styles from the active YAML config.

## Key Outputs

- `*_FiltDRR.npy`: filtered and downsampled dR/R time series, shaped as `T x H x W`.
- `*_FC_Maticx_*_excluded.npy`: ROI-level functional connectivity matrix.
- `*_correlation_matrices_*_excluded.pkl`: session-level DataFrame of FC matrices indexed by mouse and trial.
- `statistics.pkl`: summary DataFrame containing group similarity, individual similarity, relative effect, statistics, and Cohen's d.
- `Similarity/*.png`: similarity matrices, paired individual/group plots, and shuffle histograms.
- `EdgeContri/*.png` and optional matrices: circular edge contribution graphs.
- `Stat/*.png`: statistical summary plots.
- `logs/default.log`: runtime logs and statistical reports.

## Important Notes

- The code assumes Windows-style paths in the example configs, but relative paths such as `data\\processed\\` work from the repository root.
- `scripts/00_main.py` is currently empty; use the numbered scripts directly.
- `scripts/01_preprocess.py` should be run from `neuro_suite2p_env` so that `srcs/prep/suite2p.py` imports `suite2p==1.0.0.1`.
- The Suite2p wrapper has been checked against `suite2p==1.0.0.1`; registration settings such as `1Preg` and `phasecorr` are intentionally retained in the Suite2p settings block.
- Preprocessing requires a display because atlas alignment uses Matplotlib interaction.
- Suite2p is called from a temporary script using the current Python interpreter and may run for a long time on large recordings.
- `Chamber_Mask.png` is expected inside each processed trial folder before `drr_narrow()` runs. If it is not generated by your alignment workflow, place or generate it before dR/R processing.
- Numbered scripts load the active YAML config through `srcs.utils.config.load_context()`. Run scripts from the repository root so relative config paths resolve correctly.
- The project stores generated data in `data/` and logs in `logs/`; both are intentionally excluded from version control.

## Minimal Workflow Checklist

1. Create `neuro_suite2p_env` from `neuro_suite2p_env.yaml`.
2. Create `ID-mRSFC-ISOI` from `environment.yaml`.
3. Edit `configs/default.yaml`, `configs/OIS.yaml`, or create a new config for your dataset.
4. Prepare raw data in `<raw>/<dataset>/<session>/<mouse>/<trial>/`.
5. Prepare ROI mask PNG files under `paths.mask`.
6. Activate `neuro_suite2p_env`, run `01_preprocess.py`, and complete interactive atlas alignment.
7. Activate `ID-mRSFC-ISOI` and run `02_correlation.py`.
8. Confirm `Session_ids` contains the sessions you want `03_identification.py` to load.
9. Run `03_identification.py` and `04_statistics.py` in `ID-mRSFC-ISOI`.

## Citation

If this repository supports a manuscript or dataset release, add the citation information here, including article title, authors, journal/preprint server, DOI, and dataset accession if available.

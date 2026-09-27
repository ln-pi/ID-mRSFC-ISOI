# Configuration Reference

Scripts read `configs/OIS.yaml` by default. Pass another config with `-c` or `--config`:

```bash
python scripts/01_preprocess.py -c configs/default.yaml
```

## paths

```yaml
paths:
  raw: "data\\raw\\"
  processed: "data\\processed\\"
  output: "data\\output\\"
  mask: "data\\mask\\"
  log: "logs\\default.log"
```

- `raw`: root folder containing raw OIS data, read from `<raw>/<dataset>/<trial>/<session>/<mouse>/`. Each mouse folder contains `info.txt` and the recording binary files.
- `processed`: root folder for intermediate outputs.
- `output`: root folder for figures, statistics, and final outputs.
- `mask`: shared folder containing ROI PNGs and `Chamber_Mask.png` (filename case differences accepted). The chamber mask is used for every recording and excluded from ROI extraction. Mask dimensions must match the aligned images.
- `log`: runtime log path.

## variables

```yaml
variables:
  dataset: "WT"
  Session_ids: ["Session1", "Session2", "Session3", "Session4"]
  Mouse_ids: ["A01", "A02"]
  Trial_ids: ["Trial1", "Trial2"]
  Atlas: [["SM", "VIS", "AUD"], ["None"]]
```

- `dataset`: dataset/group name used under `processed` and `output`.
- `Session_ids`: sessions to process and analyze. `03_identification.py` loads pkl files in this exact order.
- `Mouse_ids`: mice to process and analyze.
- `Trial_ids`: trials to process and analyze.
- `Atlas`: ROI exclusion sets. Each inner list is converted to an atlas tag with underscores.

Examples:

- `["SM", "VIS", "AUD"]` excludes ROI names containing `SM`, `VIS`, or `AUD`.
- `["None"]` is used as the "whole brain" condition in current plotting labels.
- `["AUDv"]` excludes ROI names containing `AUDv`.

## preprocess

```yaml
preprocess:
  keepna: False
  bin:
    time: 6
    space: 4
  filter:
    highpass: 0.005
    lopass: 1
    fps: 5
  robust: 2
```

- `keepna`: whether missing/zero frames should remain as missing values. Use `False` for the current preprocessing workflow.
- `bin`: reserved for binning-related parameters.
- `filter`: intended filtering parameters. The current `Filter_2D()` implementation uses fixed band-pass values internally.
- `robust`: reserved for robust-processing settings.

## ID

```yaml
ID:
  group_mode: random
  repeat_session: null
  session_seed: 42
  shuffle: True
  n_iter: 1000
```

- `group_mode`: session grouping in `build_half_dataframe()`; any session count >= 2 is supported by the new modes.
  - `sequential`: split in `Session_ids` order. With five sessions, groups are `[1, 2, 3]` and `[3, 4, 5]`.
  - `random` (default in pipeline configs): shuffle session order once, then split equally. For odd counts, the middle session in the shuffled order is shared by both groups.
  - Legacy numeric `12`, `13`, `14` retain their original pairings for exactly four sessions only.
- `repeat_session`: optional session ID (e.g. `"Session2"`) or **one-based position** in `Session_ids`. Only valid for odd session counts in `sequential` mode. The selected session is shared; remaining sessions are split in original order. With five sessions and `repeat_session: 2`, groups are `[1, 2, 3]` and `[2, 4, 5]`. `null` shares the middle session. Direct Python calls accept the one-based position only.
- `session_seed`: nonnegative integer (default 42), shared by reproducible session grouping and permutation testing. The same config resolves to the same output directory in scripts 02–05. Null is rejected to avoid unrelated seeds between independently launched stages. `Edge.seed` independently controls within-session trial splitting.

Each half averages Fisher-z-transformed matrices equally across included sessions and trials, preserving the previous output scale. `trial_cnt=0` in direct Python calls now means all trials. The shared session contributes to both halves, so odd-session halves overlap by design.

Mouse and trial counts are not fixed at 10 and 7. All session pickle tables must contain identical mouse/trial labels in the same order and compatible ROI matrices. Change `Mouse_ids` and `Trial_ids` and regenerate correlation files when changing the dataset; identification reads the mice and trial order from those files. The averaging helper supports one or more mice/trials, but individual-vs-group comparison needs multiple mice, the current Lilliefors statistics step needs at least four mice, and edge contribution requires at least two valid trials per mouse/session and two valid mice. Statistical estimation may impose further sample-size or non-degeneracy requirements.

- `shuffle`: whether to run permutation testing in `03_identification.py`.
- `n_iter`: number of shuffle/permutation iterations.

## Edge

```yaml
Edge:
  seed: 42
  repeat: 50
  threshold: 99
  gap: 0.25
  cross_only: False
  save_mat: True
```

- `seed`: random seed for repeated edge contribution splits.
- `repeat`: number of repeated random trial splits.
- `threshold`: percentile threshold for circular graph visualization.
- `gap`: visual gap between modules in circular graphs.
- `cross_only`: whether to keep only cross-module edges.
- `save_mat`: whether to save matrix-style companion figures from circular graph outputs.

## supplement

This section controls `scripts/05_supplement.py`. The script is optional and uses the `statistics.pkl` file generated by `03_identification.py`.

```yaml
supplement:
  structure_atlases: [["SM", "VIS", "AUD"], ["PF", "LAT", "MED"], ["None"]]
  structure_group_names: ["Association", "Sensory", "Whole Brain"]
  exclude: ["AUDv"]
```

- `analysis.statistics_file`: shared statistics filename for 03, 04 and 05 (see below); the old supplement-only setting is a fallback for older configurations.
- `structure_atlases`: atlas conditions to compare in `05_supplement.py`. Each inner list is converted to the same underscore-joined row tag used in `statistics.pkl`, for example `["SM", "VIS", "AUD"]` becomes `SM_VIS_AUD`.
- `structure_group_names`: display names paired with `structure_atlases`. Keep the same length and order as `structure_atlases`.
- `exclude`: statistic rows ignored before structure-space analysis.

`05_supplement.py` extracts two metrics from each selected condition:

- `diagonal_dominance`: reported as DD.
- `normalized_distance_to_diagonal`: reported as ND.

It writes DD/ND boxplots and a DD-ND structure-space plot under `paths.output/<dataset>/Stat/`, and logs RM-ANOVA and MANOVA tables to `paths.log`.

## visualization

```yaml
visualization:
  defaults:
    dpi: 180
    show: False
    save: True
    log: False

  fc_heatmap:
    figsize: [8, 8]
    cmap: "coolwarm"
    vmin: -1
    vmax: 1

  stat_box:
    figsize: [4, 4]
    box_width: 0.1
    showfliers: False
```

Visualization presets keep plotting style out of the analysis scripts.

- `defaults`: shared render options used by scripts.
- `fc_heatmap`: ROI functional connectivity heatmap style.
- `id_shuffle_hist`: shuffle histogram style in individual identification.
- `id_paired`: paired group/individual plot style.
- `id_similarity_matrix`: similarity matrix style.
- `edge_phi`: circular graph style for group-consistency edge contribution.
- `edge_dp`: circular graph style for differential-power edge contribution.
- `stat_box`: base boxplot style.
- `stat_rel_effect_box`: extra options for relative-effect boxplots.
- `stat_group_indiv_box`: extra options for group/individual boxplots.
- `stat_mean_error`: group/individual trend-line style in statistics.
- `supplement_dd_box`: DD boxplot style in supplementary analysis.
- `supplement_nd_box`: ND boxplot style in supplementary analysis.
- `supplement_structure_space`: DD-ND structure-space scatter style.

Script-level presets are merged with plot-specific runtime values such as output names, labels, and data-dependent axis ticks.

## Choosing a Config

Use `configs/default.yaml` as a self-contained relative-path template.

Use `configs/OIS.yaml` when working with the external OIS data root currently written in that file.

For a new dataset, copy one of these files and edit:

1. `paths.raw`
2. `variables.dataset`
3. `variables.Session_ids`
4. `variables.Mouse_ids`
5. `variables.Trial_ids`
6. `variables.Atlas`


## Analysis output and resume

```yaml
analysis:
  statistics_file: "statistics.pkl"
preprocess:
  resume: True
ID:
  group_mode: random
  session_seed: 42
  repeat_session: null
Edge:
  seed: 42
```

Output layout is `<paths.output>/seed_<mode>_<session_seed>[_repeatN]_edge<edge_seed>_<configuration-hash>/<dataset>/`.
The short hash distinguishes session/data selections and ID/Edge settings. Relative paths remain relative to the launch directory. Repeating the same configuration reuses the same analysis directory. `analysis_config.json` records the configuration, actual session groups and trial order. `analysis.statistics_file` is a filename (not an absolute path); older supplement-only configurations remain supported.

03 uses the configured mouse IDs and trial order when selecting pickle entries. Each `NTrial/` directory contains identification results for the first N trials, numerical similarity/subject/statistics exports, full shuffle results, edge contributions and a one-column copy of the shared statistics pickle. A combined pickle at the analysis root serves 04/05. Edge keeps the original within-session trial split: 1Trial is explicitly unavailable; N>=2 uses only the first N trials.

04 exports ANOVA, repeated-measures ANOVA, observations and FDR-corrected post-hoc CSV tables under `Stat/Tables/`. 05 exports input observations, ANOVA tables and MANOVA tables/numerical results under `Stat/Tables/Supplement/`. Edge CSV matrices and session-wise Phi/DP pickle results are saved in addition to graph images.

With `preprocess.resume: True`, 01 validates each step's outputs before skipping it. It checks TIFF byte ranges and frame dimensions/counts, NPY structure/size, preview readability and dependencies. Step 4 skips both manual selection and transformation when its NPY and preview are complete. A running/complete JSON marker per step prevents interrupted runs being mistaken for completed runs. Changed upstream files cause dependent steps to rerun; changing the chamber mask reruns step 5. Older outputs without markers are adopted after structural validation and dependency timestamp checks. These checks detect incomplete files, not arbitrary numerical corruption; legacy output provenance cannot be reconstructed.

Set `resume: False` to rerun all five steps. Before rerunning Suite2p, an existing recording-local `suite2p/` directory is renamed to `suite2p.previous-<timestamp>/` to preserve it and avoid merging stale TIFF chunks. The legacy global `OIS_Corrected` directory is neither migrated nor deleted.

Suite2p checkpoints include the installed Suite2p version and the adapter source hash. A change to either reruns registration and downstream steps; complete legacy outputs without checkpoint metadata are still adopted after structural validation.

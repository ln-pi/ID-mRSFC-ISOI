# Configuration Reference

Use `python run.py -c configs/Aseries.yaml` or another config. With no config argument, Aseries and Nseries execute sequentially.
Relative config and data paths resolve against the P1 project root, independently of the launch directory.

## Physical recording dimensions

```yaml
variables:
  dataset: ""
  Week_ids: ["00w", "01w", "02w"]
  Mouse_ids: ["N01", "N02"]
  Run_ids: ["Run01", "Run02"]
  Atlas: [["SM", "VIS", "AUD"], ["PF", "LAT", "MED", "AUDv"], ["AUDv"]]
```

Weeks and runs always describe acquisition identities, independently of analysis.
All ID lists must be nonempty and unique. Order defines chronological week order
and cumulative run order; no lexical sorting is substituted. The former
Session_ids/Trial_ids keys must be replaced with explicit physical Week_ids/Run_ids;
old arbitrary session/trial names cannot be safely interpreted automatically.

`paths` contains raw, processed, output, mask and log roots. The dataset is appended
to raw/processed and output roots. All datasets may share `data/mask/`, including
Chamber_Mask.png and ROI PNGs. The chamber mask is excluded from ROI extraction.
Atlas lists are ROI exclusions, not inclusions.

## Input paths and preprocessing

01 and 02 search the configured week/mouse/run IDs across all six permutations
of those three directory levels, including week/run/mouse raw folders and
run/mouse/week or week/mouse/run processed folders. No data directory is renamed.
02 only accepts `<mouse>_<run>_<week>_FiltDRR.npy` (for example N01_Run01_00w_FiltDRR.npy).
01 step 5 writes the same name. Former filenames are ignored.
Multiple matching files for one recording are an error; duplicate sources are
never silently selected. Missing/conflicting/corrupt records are logged and
recorded in metadata/batch_reports while remaining records continue. Non-missing execution failures stop the workflow after that stage.

01 creates new outputs at processed/<dataset>/<week>/<mouse>/<run>/, or resumes
the single existing recording folder when present. Suite2p stays inside that
folder. Raw folders require info.txt and img*.bin. Finish copying before running
02, so incomplete files are not treated as valid input.

## Common analysis options

```yaml
analysis:
  split_axis: run          # run or week
  selection_mode: single  # single or cumulative, on the opposite axis
  statistics_file: "statistics.pkl"
ID:
  completeness_mode: strict
  group_mode: random
  repeat_id: null
  session_seed: 42
  shuffle: True
  n_iter: 1000
```

| split_axis | selection_mode | Comparison | Output conditions |
|---|---|---|---|
| run | single | Two run groups within each week | 00w, 01w, ... |
| run | cumulative | Two run groups after averaging first N weeks | 1Week, 2Week, ... |
| week | single | Two week groups using each run separately | Run01, Run02, ... |
| week | cumulative | Two week groups after averaging first N runs | 1Run, 2Run, ... |

ID first averages Fisher-z FC matrices within each half and selected opposite-axis
samples, then identifies mice. A split axis needs at least two entries.
`group_mode` is random (seeded permutation) or sequential (configured order).
Even counts split equally; odd counts share the middle entry by default.
`repeat_id` accepts an ID or one-based position on the *selected* split axis,
only for odd sequential splits. The seed retains its existing name session_seed
for continuity. Shared entries appear explicitly in the metadata.

Strict retains all configured weeks and runs, excluding any mouse with any invalid
record across that entire range and all atlases. Cover enumerates continuous
intervals in configured Week_ids order, retaining mice complete across all runs
and atlases, and maximizes mice count times week count. At least two complete mice
are required; ties prefer more mice then the earlier start. The same selected
mice and weeks are used throughout ID, shuffle, Edge and statistics. If a selected
week split has fewer than two weeks, 03 stops rather than inventing another group.
Missing whole week/atlas PKLs or incomplete latest 02 reports stop 03 explicitly.

Identification/cohort_selection.json records candidates, scores, exclusions and
all actual condition selections. The aggregate and numeric results carry week,
run, mouse, split-group and selection-mode metadata. 04/05 check this metadata
against the current config before processing. 04 labels its repeated factor week
or run and uses real mouse IDs; a single condition uses the condition factor only.

## Edge

Edge shares analysis.split_axis and the actual ID groups; there is no Edge.split_axis.
It computes raw-FC Phi/DP separately within each selected opposite-axis entry, then
averages those values for the figure. Thus run/single can compute Edge from two
runs in one week, and week/single can compute Edge from multiple weeks of one run.
At least two split-axis samples and two complete mice are required.

The shared groups are fixed. Edge.repeat does not repeat identical pairs, and
Edge.seed does not determine an independent split in this workflow. These options
remain for compatibility with direct legacy Edge calls and participate in config
isolation. This change in pairing means old Edge outputs are not expected to be
numerically reproduced. ID with run splits preserves the former run-as-session
identification computation, apart from explicit output/metadata labels.

Edge.threshold, gap, cross_only and save_mat control figures as before.

## Outputs and configuration isolation

The complete parsed config determines a stable short hash, with output root:
`<output>/seed_<group_mode>_<session_seed>[_repeatID]_edge<edge_seed>_<hash>/<dataset>/`.
Changing any config value, including analysis axes, selects a different directory.
Set the final config before running 02–05. No cross-config cache fallback or
automatic old-PKL migration is performed. Source FC computations are invariant
under input directory permutation, but table schema and config identity changed.

```text
Correlation/
  FC_Matrix/<week>_<mouse>_<run>_FC_Matrix_<atlas>_excluded.{npy,cache.pkl}
  Figures/
  <week>_correlation_matrices_<atlas>_excluded.pkl
Identification/
  cohort_selection.json
  statistics.pkl
  <condition>/Numeric/
  <condition>/EdgeContri/
  <condition>/Similarity/
  <condition>/Stat/
Statistics/
  analysis_metadata.json
  Tables/
    Supplement/
BatchReports/
analysis_config.json
```

02 PKLs always have mouse rows and run columns, one table per week/atlas. Missing
cells remain None. Reports use week, mouse, run, atlas identities and include
actual source paths for successful recordings. FC reuse validates source metadata,
masks, computation signature and config identity; missing sources cannot be hidden
by a cache. Plot settings/readability are checked independently of FC reuse.

04 writes observations, ANOVA and FDR post-hoc tables. 05 uses
supplement.structure_atlases, structure_group_names and exclude, exporting DD/ND
observations with actual analysis-condition labels. The existing 05 inferential
model still treats condition-level summaries as repeated observations; the pending
statistical-design concerns described in docs/optimization.md are not resolved by
this dimension migration. Conditions are not automatically independent samples. With only one condition,
05 saves descriptions/plots and inference_status.json with unavailable status.

## Preprocessing resume

With `preprocess.resume: True`, 01 validates each step's outputs before skipping it. It checks TIFF byte ranges and frame dimensions/counts, NPY structure/size, preview readability and dependencies. Step 4 skips both manual selection and transformation when its NPY and preview are complete. A running/complete JSON marker per step prevents interrupted runs being mistaken for completed runs. Changed upstream files cause dependent steps to rerun; changing the chamber mask reruns step 5. Older outputs without markers are adopted after structural validation and dependency timestamp checks. These checks detect incomplete files, not arbitrary numerical corruption; legacy output provenance cannot be reconstructed.

Set `resume: False` to rerun all five steps. Before rerunning Suite2p, an existing recording-local `suite2p/` directory is renamed to `suite2p.previous-<timestamp>/` to preserve it and avoid merging stale TIFF chunks. The legacy global `OIS_Corrected` directory is neither migrated nor deleted.

Suite2p checkpoints include the installed Suite2p version and the adapter source hash. A change to either reruns registration and downstream steps; complete legacy outputs without checkpoint metadata are still adopted after structural validation.


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


## P2 FC reuse and output layout

Use `--start identification --import-fc SOURCE_OUTPUT` with `-c` to import P2 FC tables. Source output root or Correlation folder is accepted; raw/processed inputs are not required for imported FC analysis. The complete parsed config still determines analysis_tag. Relative path resolution does not rewrite the config hash. Output categories are Correlation, results/combined, report, and metadata; see [workflow](workflow.md).


## Configuration templates and similarity color scales

`Aseries.yaml` and `Nseries.yaml` are runnable configurations. `default.yaml` retains the same property tree with all leaf values empty (YAML null); fill it before use. The former OIS.yaml has been removed. Legacy numbered scripts default to Aseries.yaml.

`visualization.id_similarity_matrix.limits_by_atlas` sets fixed Fisher-z similarity limits shared across A/N and Run01/Run02: SM_VIS_AUD `[1.3, 2.9]`, PF_LAT_MED_AUDv `[0.9, 2.7]`, AUDv `[1.0, 2.8]`. Raw FC heatmaps retain `[-1, 1]`. Full-config hashing means these visualization changes produce new output tags. Existing figures are not rewritten automatically.

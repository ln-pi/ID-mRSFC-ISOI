# ID-mRSFC-ISOI

A Python workflow for studying individual differences and longitudinal changes in mouse resting-state functional connectivity (mRSFC) using intrinsic signal optical imaging (ISOI/OIS).

The project connects acquisition-level preprocessing, cortical functional connectivity (FC), individual identification, and longitudinal network analysis. Each mouse, week, and recording Run remains identifiable throughout the workflow, supporting both comparisons of individual recordings and comparisons of connectivity patterns averaged across observations.

## Features

- OIS binary recording conversion, Suite2p motion correction, and interactive atlas registration.
- Reflectance-change preprocessing, detrending, global signal regression, filtering, and downsampling.
- ROI time-series extraction and labeled cortical FC matrices.
- Shared FC edge arrays for combined and longitudinal analyses.
- Within-Run cross-week similarity, phase-based comparisons, and individual identifiability metrics.
- Euclidean/RMS and correlation distances, group dispersion, and individual deviation examples.
- Leave-one-out Delta FC maps, edge variance, regional/network summaries, and deviation direction.
- Same-day Run repeatability and changes from each mouse's reference week.
- Edgewise stable-mouse ICC, absolute variance components, and stable individual effect maps.
- Delta FC persistence, animal-block trajectory intervals, and network random slopes.
- Identification permutation analysis, edge contribution, and statistical summaries.
- Figures generated independently from saved analysis results.

## Data organization

The recording structure is **mouse × week × Run**. The example A-series and N-series configurations cover `00w–13w` with two same-day recordings, `Run01` and `Run02`. Mouse lists, paths, cortical partitions, and analysis selections are specified in YAML configuration files.

Each cortical partition has its own edge arrays:

```text
F_r[mouse, week, run, edge]    Original Pearson FC edge weights
F_z[mouse, week, run, edge]    Fisher-z-transformed FC edge weights
valid_record[mouse, week, run]
```

The edge dimension contains each undirected, off-diagonal ROI connection once, in a fixed strict-upper-triangle order. Mouse, week, Run, ROI, and edge labels accompany the arrays. Missing recordings retain their positions and are represented by masks and NaN values.

The example cortical partitions are configured through ROI exclusion lists:

| Configuration label | Retained cortical partition |
| --- | --- |
| `SM_VIS_AUD` | Association cortex |
| `PF_LAT_MED_AUDv` | Sensory cortex, excluding AUDv |
| `AUDv` | Whole cortex, excluding AUDv |

## Analysis methods

### Functional connectivity

The preprocessing workflow converts raw recordings into motion-corrected, atlas-aligned reflectance-change signals. ROI time series are extracted for each recording separately, and Pearson correlations between ROI signals form the FC matrix. FC edges are retained on both the original correlation scale and the Fisher-z scale.

### Longitudinal pairwise analysis

For each cortical partition and each Run, Pearson correlation across Fisher-z FC edges defines a four-dimensional pattern-similarity array:

```text
S[mouse_i, mouse_j, week_u, week_v]
```

`S[:, :, u, v]` compares all mice at week u with all mice at week v. `S[i, i, :, :]` describes the cross-week stability of one mouse's connectivity pattern. Reverse week-pair panels are transposes of one another; individual cross-week panels need not be symmetric.

Pattern similarity is saved as Pearson r and Fisher-z values. Comparisons of a recording with itself are excluded from statistical summaries. Missing recordings and constant edge patterns produce undefined similarities rather than zero values.

Time windows select observations from the full arrays. The example configuration includes onset (`00w–02w`), a reference week (`03w`), early (`04w–06w`), middle (`07w–10w`), and late (`11w–13w`) windows. N-series windows carry the study's disease-time interpretation; A-series windows are aligned observation windows. `03w` is a post-injection reference. Within-window analyses use distinct week pairs, and between-window analyses retain all selected directed week pairs.

### Combined identification

Combined analysis divides the selected observation axis into two non-overlapping groups, averages Fisher-z FC edges within each group, and identifies mice by comparing the resulting connectivity patterns. The example setup splits 14 weeks into two random groups of seven, separately for each Run.

This calculation averages FC edges **before** computing pattern similarity. It estimates a different quantity from averaging cross-week similarities. Combined similarity is represented on a Fisher-z scale capped at ±4. Group selection, random seed, and accumulation settings are configurable.

For each comparison, the project reports:

- Own-match similarity and mean nonself similarity.
- Self-match advantage: own-match similarity minus mean nonself similarity.
- Row and column matching percentiles, with ties contributing half a match.
- Identification accuracy and its permutation distribution.
- Edge contributions: group consistency (Phi) and differential power (DP).

### Repeatability and variability

Descriptive analysis measures FC differences in Fisher-z edge space. For E edges, Euclidean distance is `sqrt(sum((f_i-f_j)^2))` and RMS distance is Euclidean distance divided by `sqrt(E)`. Pattern distance is `1-Pearson(f_i,f_j)`; it is reported separately from amplitude distances.

Each mouse's Delta FC is its edge vector minus the mean of the other mice in the same group, week, and Run view. Its RMS summarizes the magnitude of deviation. Group dispersion uses the full group mean and denominator N; edgewise sample variance uses N−1. Positive Delta FC means stronger connectivity relative to the reference, and negative Delta FC means weaker connectivity. These signs do not indicate excitation or inhibition.

The default views include each Run separately and FC averaged across all configured Runs **before** calculating metrics. The averaged view requires every Run to be valid. Results include all mice available each week and a cohort complete across the selected window. The largest LOO RMS identifies a descriptive example; it does not label or exclude an animal as a statistical outlier.

Delta FC and variance are summarized by ROI, network pair, and hemisphere pair. ROI contributions split each edge's squared deviation or variance equally between its two endpoints; network pairs retain unique edges. Both normalized contributions and per-edge averages are saved. The built-in network grouping follows the atlas's anatomical families; explicit ROI-to-network mappings can define other networks.

Same-day outputs quantify pattern correlation, RMS differences, and signed shifts between Runs. Changes from each mouse's own reference week remain separate from same-week group Delta FC. Amplitude/topology outputs remain descriptive; variance models and trajectory uncertainty are separate stages.

### Stability and longitudinal evolution

Stability analysis compares raw FC and Delta FC between same-day Runs using common reference animals. Across weeks, each mouse's Delta FC vectors are correlated over edges; the same-record diagonal is excluded. Dynamic and fixed-window references are reported separately.

For each edge, a balanced repeated-measures model includes fixed week and Run effects, a stable mouse effect, a mouse-by-week effect, and Run-level residual error. Stable-mouse ICC is `mouse_variance / (mouse_variance + mouse_week_variance + error_variance)`. Absolute variance, animal-bootstrap intervals, boundary status, and conditional stable-mouse effect maps accompany ICC. This model uses animals complete across the selected weeks and Runs, with an explicit cohort audit. The stable effect may include persistent imaging or anatomical differences.

Longitudinal analysis reports group dispersion and individual deviation trajectories, rank persistence, and dispersion trend intervals obtained by resampling whole animals. Available and fixed cohorts remain separate; increasing/decreasing dispersion labels describe an overall time trend rather than a causal or necessarily monotonic process.

Network FC means are also fitted with linear time and Run fixed effects, correlated mouse intercepts/slopes, and a mouse-by-week random effect using [Statsmodels MixedLM](https://www.statsmodels.org/stable/generated/statsmodels.regression.mixed_linear_model.MixedLM.html). Outputs include conditional individual slopes, variance estimates, convergence/boundary flags, optimizer warnings, and residual diagnostics. Nonconverged fits do not produce individual slope predictions. Informative loss of severely affected animals is not resolved by these models.

Pairwise summaries provide both mice available at each week pair and a fixed cohort complete across the selected window. Candidate sets are recorded, and matching metrics are computed within each set. Combined analyses use a separately configured completeness policy.

Statistical summaries include paired effect sizes, repeated-measures analysis, and FDR-adjusted post-hoc comparisons. Mice are the experimental units; weeks and Runs are repeated observations. Phase-level pairwise and variability outputs describe connectivity patterns and trajectories.

## Workflow

| Stage | Purpose |
| --- | --- |
| `preprocess` | Convert raw recordings, perform motion correction and atlas alignment, and prepare filtered signals |
| `fc` | Compute labeled FC matrices for individual recordings |
| `edges` | Assemble the shared mouse × week × Run × edge arrays |
| `pairwise` | Calculate per-Run cross-week similarities and window-based mouse metrics |
| `descriptive` | Compute distances, dispersion, Delta FC, edge/regional variability, and Run repeatability |
| `stability` | Estimate Run repeatability, Delta FC persistence, stable components, and edgewise ICC |
| `longitudinal` | Calculate trajectory intervals, dispersion trends, rank persistence, and network random slopes |
| `combined` | Calculate identification, permutation distributions, and edge contributions from grouped observations |
| `statistics` | Export combined-analysis subject metrics and statistical tables |
| `report` | Generate figures for saved analysis sections |

## Installation and configuration

Create and activate the supplied Conda environment:

```powershell
conda env create -f environment.yaml
conda activate ID-mRSFC-ISOI
```

Before running, configure the data and mask paths in `configs/Aseries.yaml`, `configs/Nseries.yaml`, or `configs/default.yaml`. Relative paths resolve from the project root. Preprocessing requires the raw recording files and atlas-registration interaction; downstream analyses can also start from existing labeled FC results.

See the [configuration reference](docs/config.md) for input conventions, analysis parameters, and array contents, and the [workflow guide](docs/workflow.md) for acquisition layout and preprocessing details.

## Running the project

Run commands from the project root in PowerShell:

```powershell
# Inspect the execution plan and available inputs.
python run.py --dry-run

# Run A-series and N-series from FC computation onward.
python run.py --configs configs/Aseries.yaml configs/Nseries.yaml

# Run preprocessing for one configuration.
python run.py -c configs/Aseries.yaml --stage preprocess
```

`--stage` runs one stage. `--start` runs that stage and all subsequent stages in the order listed above. The default starts at `fc`, processing A-series followed by N-series. The PowerShell wrapper accepts the same arguments, for example `.\run.ps1 --start edges`.

To use an existing FC output directory:

```powershell
python run.py -c configs/Aseries.yaml --start edges --import-fc "F:\path\to\saved_FC_output"
```

The source can be an output root or its `correlation` directory. Recording identities, ROI labels, source information, and missing observations accompany the imported matrices.

To run individual analyses or redraw saved results:

```powershell
python run.py -c configs/Aseries.yaml --stage pairwise
python run.py -c configs/Aseries.yaml --stage descriptive
python run.py -c configs/Aseries.yaml --stage stability
python run.py -c configs/Aseries.yaml --stage longitudinal
python run.py -c configs/Aseries.yaml --stage combined
python run.py -c configs/Aseries.yaml --stage report
python run.py -c configs/Aseries.yaml --stage report --report-section combined
python run.py -c configs/Aseries.yaml --stage report --report-section descriptive
python run.py -c configs/Aseries.yaml --stage report --report-section stability
python run.py -c configs/Aseries.yaml --stage report --report-section longitudinal
```

`--report-section` selects `pairwise`, `descriptive`, `stability`, `longitudinal`, `combined`, `statistics`, or `all`. Numerical outputs and figure generation are separate, so saved analyses can be redrawn without rerunning their calculations.

## Outputs

Results are organized by configuration:

```text
data/output/<analysis_tag>/<dataset>/
├── correlation/             Labeled FC matrices
├── results/
│   ├── edges/               FC edge arrays, labels, and recording masks
│   ├── pairwise/            Per-Run similarity arrays and week-pair/window summaries
│   ├── descriptive/v2/      Distances, dispersion, Delta FC, regions, and examples
│   ├── stability/           Repeatability, Delta FC persistence, variance components, ICC
│   ├── longitudinal/        Trajectories, uncertainty, trends, network slopes and diagnostics
│   ├── combined/            Identification, permutations, and Phi/DP
│   └── statistics/          Subject metrics and statistical tables
├── report/                  Figures and report indexes
└── metadata/                Configuration, source information, and recording/stage records
```

Logs are stored at the configured log paths. Batch execution plans and status summaries are written under `reports/automatic/`. Generated data, reports, and logs are local outputs.

## Project structure

```text
configs/                     Example study configurations
docs/                        Workflow and configuration documentation
srcs/
├── analysis/                Connectivity, identification, and longitudinal methods
├── preprocessing/           Recording readers, signal processing, and registration
├── visualization/           Diagnostic figures and saved-result rendering
├── workflow/                Analysis stages and report dispatch
├── utils/                   Configuration, data access, output handling, and provenance
└── atlas/                   Cortical atlas resources
run.py                       Python command-line entry point
run.ps1                      PowerShell entry point
environment.yaml             Conda environment
pyproject.toml               Package metadata and dependencies
```

## Data availability

Raw experimental recordings are not distributed with this repository. Contact the authors regarding data access and collaboration.

## Citation and acknowledgments

Please cite the software using [CITATION.cff](CITATION.cff). Related literature includes:

- Finn, E. S., et al. (2015). Functional connectome fingerprinting: identifying individuals using patterns of brain connectivity. *Nature Neuroscience*. [DOI: 10.1038/nn.4135](https://doi.org/10.1038/nn.4135)
- White, B. R., et al. (2011). Imaging of functional connectivity in the mouse brain. *PLOS ONE*. [DOI: 10.1371/journal.pone.0016322](https://doi.org/10.1371/journal.pone.0016322)
- Bergmann, E., et al. (2020). Individual variability in functional connectivity architecture of the mouse brain. *Communications Biology*. [DOI: 10.1038/s42003-020-01472-5](https://doi.org/10.1038/s42003-020-01472-5)

Circular connectivity visualization is adapted from [paul-kassebaum-mathworks/circularGraph](https://github.com/paul-kassebaum-mathworks/circularGraph).

## License

This project is distributed under the [MIT License](LICENSE).

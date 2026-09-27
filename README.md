# ID-mRSFC-ISOI

This repository contains a Python workflow for mouse resting-state functional connectivity analysis from intrinsic signal optical imaging (ISOI/OIS) data. The pipeline converts raw OIS binary recordings to TIFF stacks, performs Suite2p motion correction and atlas registration, extracts ROI-level time series, computes functional connectivity (FC) matrices, and evaluates individual-level network identifiability from mRSFC patterns.

The code is organized as research scripts rather than a packaged command-line application. Most runtime behavior is controlled by YAML configuration files in `configs/`.

## Documentation

- Detailed workflow: [docs/workflow.md](docs/workflow.md)
- Configuration reference: [docs/config.md](docs/config.md)
- Citation metadata: [CITATION.cff](CITATION.cff)
- License: [LICENSE](LICENSE)

## Overview

Main analysis stages:

1. Preprocess raw OIS data with `scripts/01_preprocess.py`.
2. Compute ROI-level FC matrices with `scripts/02_correlation.py`.
3. Run individual identification and edge-contribution analysis with `scripts/03_identification.py`.
4. Run statistical summaries with `scripts/04_statistics.py`.
5. Optionally run supplementary DD/ND structure analysis with `scripts/05_supplement.py`.

Core capabilities:

- OIS200-style raw binary reading.
- TIFF conversion and Suite2p motion correction.
- Interactive atlas alignment.
- dR/R, detrending, global signal regression, filtering, and downsampling.
- ROI time-series extraction and FC matrix estimation.
- Individual identification, permutation testing, edge contribution, and statistical plotting.

Repository layout:

```text
ID-mRSFC-ISOI/
|-- configs/                 # YAML runtime configs
|-- data/                    # Local raw, mask, processed, and output data
|-- docs/                    # Workflow and config documentation
|-- logs/                    # Runtime logs
|-- notebooks/               # Exploratory notebooks
|-- scripts/                 # Numbered pipeline scripts
|-- srcs/                    # Core Python modules
|-- tests/environment.yaml   # Unified Python 3.13.13 analysis + Suite2p environment
|-- CITATION.cff
|-- LICENSE
`-- README.md
```

Generated data and logs are not version-controlled.

## Quick Start

Create the unified environment:

```bash
conda env create -f tests/environment.yaml
```

Prepare your active config before running. `configs/default.yaml` is the recommended relative-path template; `configs/OIS.yaml` contains a lab-machine raw-data path example. See [docs/config.md](docs/config.md) for all fields.

Run preprocessing and analysis in the same environment:

```bash
conda activate ID-mRSFC-ISOI_suite2p
python scripts/01_preprocess.py -c configs/default.yaml
```

Run downstream analysis in the same environment:

```bash
conda activate ID-mRSFC-ISOI_suite2p
python scripts/02_correlation.py -c configs/default.yaml
python scripts/03_identification.py -c configs/default.yaml
python scripts/04_statistics.py -c configs/default.yaml
```

Optional supplementary analysis:

```bash
python scripts/05_supplement.py -c configs/default.yaml
```

For required input layout, masks, interactive atlas alignment, expected outputs, and practical notes, see [docs/workflow.md](docs/workflow.md).

## Known Issues

Some implementation details are intentionally left unchanged for now:

- `keepna=True` in raw binary frame reading may not preserve missing frames as NaN because adjusted frames are stored in an unsigned integer array.
- `srcs/atlas/Atlas_Mask.py::Get_Func_Mask()` appears to contain a non-main-path indexing issue and should be checked before using that helper directly.

See [docs/workflow.md](docs/workflow.md) for runtime prerequisites and [docs/config.md](docs/config.md) for configuration details.

## Data Availability

Raw experimental data and lab-internal core code are not publicly distributed with this repository. Please contact the authors for data access or collaboration requests.

## Citation

If you use this repository, cite the software metadata in [CITATION.cff](CITATION.cff) and the core references below as appropriate.

Core literature:

- Finn, E. S., et al. (2015). Functional connectome fingerprinting: identifying individuals using patterns of brain connectivity. Nature Neuroscience. DOI: [10.1038/nn.4135](https://doi.org/10.1038/nn.4135)
- White, B. R., et al. (2011). Imaging of functional connectivity in the mouse brain. PLOS ONE. DOI: [10.1371/journal.pone.0016322](https://doi.org/10.1371/journal.pone.0016322)
- Bergmann, E., et al. (2020). Individual variability in functional connectivity architecture of the mouse brain. Communications Biology. DOI: [10.1038/s42003-020-01472-5](https://doi.org/10.1038/s42003-020-01472-5)

Code references:

- `circularGraph` visualization is adapted from [paul-kassebaum-mathworks/circularGraph](https://github.com/paul-kassebaum-mathworks/circularGraph).
- Additional lab core code is internal and not publicly released.

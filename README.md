# Analysis
[![Build Status](https://dev.azure.com/pshama/analysis/_apis/build/status/uppittu11.analysis?branchName=master)](https://dev.azure.com/pshama/analysis/_build/latest?definitionId=1&branchName=master)
[![codecov](https://codecov.io/gh/uppittu11/analysis/branch/master/graph/badge.svg)](https://codecov.io/gh/uppittu11/analysis)

Analysis scripts for molecular dynamics simulations of stratum corneum lipid multilayers.

## Installation Instructions
The package can be installed via `pip`

1. Clone the repository
```bash
git clone https://github.com/uppittu11/analysis.git
```

2. Install the package
```bash
cd analysis
pip install -e .
```

3. Run tests for the package
```bash
pytest -v analysis
```

I suggest you install this package in a separate `conda` environment so the `analyze` command does not clobber anything in your default namespace

## Usage
To calculate basic structural properties (tilt angle, repeat distance, nematic order, and area per lipid) using a parallel multiprocessing pool, use the bash command `analyze`.

### Command line arguments for `analyze`
`-f` trajectory file (needs to be a format loadable by MDTraj)

`-c` topology file (needs to be a format loadable by MDTraj)

`-o` output directory

`-n` number of leaflets

`--cg` use this flag to denote a coarse-grain system

`--reload` use this flag to ignore cached trajectory (typically located in the output directory)

`--min` minimum z position of lipid headgroups to include in analysis

`--max` maximum z position of lipid headgroups to include in analysis. Use `--min` and `--max` to specify one or more layers in a multilayer system

`--center` use this flag to center the bilayer at the origin


### Loading in results
`analysis` will output a `pickle`d dictionary file containing structural properties calculated for each frame of the trajectory (default: `outputDirectory/results.p`. This output file can be loaded and analyzed using functions in the `analysis.data` module.

## Currently Supported Measurements
- Tilt Angles
- Bilayer/Multilayer Height
- Nematic Order Parameter
- Area per Lipid

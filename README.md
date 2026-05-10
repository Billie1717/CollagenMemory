# Basement membrane sets the timescale of tissue mechanical memory — simulation code

This repository contains the molecular dynamics code, analysis scripts and plotting
notebooks used to produce the simulation results in:

> Khalilgharibi, N., Meadowcroft, B., Šarić, A., Mao, Y. (2026).
> *Basement membrane sets the timescale of tissue mechanical memory.*
> bioRxiv. https://doi.org/10.64898/2026.04.23.720416

The code implements a coarse-grained model of the collagen IV network of the
basement membrane (BM), simulates uniaxial stretch and annular ablation
experiments, and reproduces the alignment, stress and connectivity analyses
shown in Figures 3, 4 and S3.

---

(Billie to do):
1. Make quickstart usable
2. Better instructions for 'running the simulations'
3. Better instructions for Data processing

## Table of contents

1. [Overview of the model](#overview-of-the-model)
2. [Repository layout](#repository-layout)
3. [Requirements](#requirements)
4. [Quick start](#quick-start)
5. [Running the simulations](#running-the-simulations)
6. [Data processing](#data-processing)
7. [Reproducing the figures](#reproducing-the-figures)
8. [Data conventions](#data-conventions)
9. [Citing this work](#citing-this-work)
10. [Contact](#contact)

---

## Overview of the model

Each collagen IV protomer is represented as a two-bead rod: one bead for the
NC1 domain and one for the 7S domain. NC1 beads bind pairwise; 7S beads can
bind up to three other 7S beads, so protomers self-assemble into a percolating
network resembling the BM. Bond formation and breaking is implemented in
LAMMPS via the `bond/react` package, which requires one `pre` / `post` / `map`
template triple per allowed change of NC1 / 7S valency — hence the large
`templates/` directories described below.

Key simulation stages are:

1. **Self-assembly + equilibration** of the network from a random initial
   configuration (`Initial_config/`).
2. **Uniaxial stretch** of the equilibrated network at a chosen strain.
3. (Optional) **Bond rearrangement** — stochastic bond breaking and reforming
   without adding new protomers — to model network remodelling.
4. **Annular ablation** in silico, used to calibrate simulation timescales
   against the experimental ablation recoil.
5. **Mutant networks**: a fraction of 7S beads has reduced valency, mimicking
   the 7S mutation. This is implemented by introducing a new bead type
   (type 9) which represents 7S ends that cannot bind. Type 9 does not
   participate in any binding reactions, in contrast to type 2 (the
   "normal" 7S binding ends that have not yet bound).

Outputs from the simulations are time series of stress, molecular alignment
(`MolAlign*`), order parameter snapshots (AFT analogue) and LAMMPS
trajectories.

A simplified, heavily-commented, **non-runnable** companion script
`model_scripts/Stretch/build_Stretch_SimpleAnnotated.py` walks through the
real `build_Stretch.py` block by block and is the recommended starting point
for understanding what the LAMMPS input file actually does.

> TODO: short reference to the specific equations / potentials, plus the
> previous model paper this builds on (ref. 37 in the manuscript).

---

## Repository layout

```
GithubCode/
├── README.md
├── MyLooseInstructions.md            # Author's working notes (cluster paths, TODOs)
├── lammps_src/                       # Custom LAMMPS REACTION package source
│   ├── fix_bond_react.cpp
│   └── fix_bond_react.h
├── model_scripts/
│   ├── Equilibrate.sh                # Top-level wrapper: equilibrate a network
│   ├── Initial_config/               # Generate initial network configs
│   │   ├── Readme.md
│   │   ├── create_dumbbell_lattice_slabNC1_type9.py
│   │   ├── create_7s_configs.py
│   │   ├── dumbbells_lattice_slab_n15680_bl3.00_Lx164.00_rho5.000e-02.lammpsdata
│   │   └── dumbbells_7s_frac0.44_n15680_bl3.00_Lx164.00_rho5.000e-02.lammpsdata
│   ├── Stretch/                      # Fig. 3 — WT stretch ± remodelling
│   │   ├── build_Equilibrate.py
│   │   ├── build_Stretch.py                  # the working script
│   │   ├── build_Stretch_SimpleAnnotated.py  # commented walk-through (read me first)
│   │   ├── build_StretchNoRemodel.py
│   │   ├── build_StretchRESTART.py           # restart a relaxation after stretch
│   │   ├── build_StretchNoRemodelRESTART.py
│   │   ├── Stretch.sh                        # SLURM wrapper for build_Stretch.py
│   │   ├── Stretch_NoRemodel.sh
│   │   ├── StretchRESTART.sh
│   │   ├── Stretch_NoRemodelRESTART.sh
│   │   └── templates/                        # bond/react templates (pre/post/map)
│   │       ├── NC1_6N/                       # NC1 valency-change reactions
│   │       └── OpenTetramerToRhomboidTetramer/
│   ├── 7s_Perturbations/             # Fig. 7 — 7S-mutant networks
│   │   ├── build_Equilibrate.py
│   │   ├── build_Stretch.py
│   │   ├── Equilibrate.sh
│   │   └── Stretch.sh
│   ├── Ablations/                    # Fig. S3 — annular ablations on large networks
│   │   ├── build_StretchAblate.py
│   │   ├── Stretch_Ablation.sh
│   │   ├── Stretch_NoAblation.sh             # benchmarking, Fig. S3a
│   │   └── SimulationsStartFromRestartFile/
│   │       └── run_data4_T1.465000.restart   # equilibrated large network used as input
│   ├── GaussianSpread/               # Render LAMMPS frames as AFT-analysable images
│   │   ├── readme.md
│   │   ├── PointSpreadSingle.py              # WT
│   │   ├── PointSpreadSingle_7s.py           # 7S mutant
│   │   ├── MetaScript.py
│   │   ├── MetaScript_7s.py
│   │   ├── ImageCreation.py
│   │   ├── CreateImages.sh
│   │   ├── CreateImages_7s.sh
│   │   ├── GaussianSpreadImages.sh
│   │   └── GaussianSpreadImages_7s.sh
│   └── DataProcessing/               # Stress / alignment / type-counting analyses
│       ├── MoleculeAngle.py
│       ├── MoleculeAlignment.sh
│       ├── Histogram_Types.py
│       └── HistogrammingTypes.sh
└── Figures/                          # One folder per figure in the paper
    ├── DataFrom.txt                  # Provenance: cluster paths the data came from
    ├── FigAlignment/                 # Fig. 3 — alignment & stress, WT network
    │   ├── Alignment.ipynb
    │   ├── Mechanics.ipynb
    │   ├── MergingFiles.ipynb        # Concatenates restart segments (maintenance)
    │   ├── Data/
    │   │   ├── Equilibration/
    │   │   ├── NoStretch/
    │   │   ├── NoRemodel_Xstretch100/
    │   │   └── Remodel_Xstretch100/
    │   └── Plots/
    ├── FigAblations/                 # Fig. S3 / calibration — simulated annular cuts
    │   ├── PlottingAblations.ipynb
    │   └── Data_all/
    ├── Fig7sMutants/                 # Fig. 7 — reduced 7S valency networks
    │   ├── Alignment.ipynb
    │   ├── Mechanics7S.ipynb
    │   └── Data/
    │       ├── Equilibration_Fraction0.56/
    │       ├── Equilibration_Fraction0.75/
    │       ├── XStretch100_Fraction0.56/
    │       └── XStretch100_Fraction0.75/
    └── Snapshots/                    # Rendered network snapshots (OVITO scenes + tif/png)
        ├── FilesFrom.md
        ├── Simulation_Remodel_All.ovito
        ├── Fig3_Simulation_Box.png
        ├── Fig3GaussianSnapshots/
        └── Fig7sMutantGaussian/
```

---

## Requirements

**Simulations**

- **LAMMPS** version `lammps-15Jun2023`. LAMMPS must be **built with the
  custom `REACTION` package source files in `lammps_src/`** — replace
  `src/REACTION/fix_bond_react.{cpp,h}` in the LAMMPS source tree with the
  versions in this repo, then build as usual.
- A C++ compiler (for building LAMMPS).
- MPI (optional) for parallel runs.

**Analysis & plotting**

- Python ≥ 3.9
- `numpy`, `scipy`, `pandas`, `matplotlib`, `seaborn`, `jupyter`
- AFT (Alignment by Fourier Transform) workflow — Marcotti, S. et al. (2021) 10.3389/fcomp.2021.745831
- OVITO (for opening the `.ovito` snapshot scenes)

> TODO: add an `environment.yml` or `requirements.txt`.

---

## Quick start

```bash
# 1. clone
git clone <repo-url>
cd GithubCode

# 2. build LAMMPS 15Jun2023 with the custom REACTION src in lammps_src/
#    (see lammps_src/ — replace src/REACTION/fix_bond_react.{cpp,h} before building)

# 3. generate or pick an initial network configuration
#    (see model_scripts/Initial_config/Readme.md for prompted parameters)

# 4. equilibrate the network
cd model_scripts
bash Equilibrate.sh
# under the hood:
#   mpirun -np 1 /path/to/lammps-15Jun2023/src/lmp_mpi -in collagen.in
# the networks in the paper were equilibrated for 3e6 timesteps

# 5. stretch + relax (with bond rearrangement = remodelling)
cd Stretch
bash Stretch.sh

# 6. analyse and plot
cd ../../Figures/FigAlignment
jupyter lab Alignment.ipynb
```

Each notebook expects to be run from its own folder, with the `Data/` folder
sitting next to it. The data paths inside the notebooks are currently
hard-coded as
`/Users/billiemeadowcroft/Dropbox/Collagen/NargessPlotting/GithubCode/...`
— **TODO: replace these with relative paths so the notebooks are portable.**

---

## Running the simulations

Each `model_scripts/<Stage>/` folder contains the LAMMPS inputs and SLURM
wrappers for one experiment. The general pattern is: a top-level `*.sh`
launcher reads parameters, calls a `build_*.py` script that writes a
LAMMPS input deck (`collagen.in`) plus a SLURM submission script
(`runscript.sh`), and then submits the job.

### Stretch (Fig. 3)

The main "stretch + remodelling" experiment.

1. Equilibrate a network with `Equilibrate.sh` — change the path to point
   to one of the initial configurations in `Initial_config/dumbbells*`.
2. Run equilibration:
   ```bash
   mpirun -np 1 /path/to/lammps-15Jun2023/src/lmp_mpi -in collagen.in
   ```
   The networks in the paper were equilibrated for **3e6** timesteps.
3. Stretch and relax, with or without remodelling:
   - With remodelling: `Stretch.sh` → `build_Stretch.py`
   - Without remodelling: `Stretch_NoRemodel.sh` → `build_StretchNoRemodel.py`
4. The `*RESTART.sh` / `build_*RESTART.py` variants restart a relaxation from
   a LAMMPS restart file. In the paper this was used multiple times back-to-
   back to get long enough trajectories to see networks relax their alignment
   via remodelling.

The output of step 3/4 is then ready for data processing with `GaussianSpread`
(see below).

To **understand** the stretch + remodelling code, read
`build_Stretch_SimpleAnnotated.py` first. It is a commented, simplified
version of `build_Stretch.py`, intended for reading rather than running
(it is not a working script).

### 7S perturbations (Fig. 7)

The 7S-perturbation simulations have different molecule templates and
molecule types from the WT case. A new bead type (type 9) is introduced
to represent 7S ends that cannot bind: type 9 is excluded from all binding
reactions, in contrast to type 2 (the "normal", as-yet-unbound 7S ends).
A different initial configuration must therefore be used — one in which a
chosen fraction of molecules carry these "perturbed" type-9 ends. Generate
these with `model_scripts/Initial_config/create_7s_configs.py`; the
templates live in `model_scripts/Stretch/templates/`.

### Ablations (Fig. S3)

Ablation simulations are run on much larger networks. The equilibrated input
is provided as a LAMMPS restart file in
`model_scripts/Ablations/SimulationsStartFromRestartFile/`. The networks are
stretched, relaxed for a short time, and then ablated (`Stretch_Ablation.sh`).
For benchmarking the stress in the absence of an ablation (Fig. S3a),
`Stretch_NoAblation.sh` runs the same protocol without the ablation step.
The bond/react templates are the same as those used for the Stretch
experiments.

### Understanding the bond/react templates

Bond formation and breaking of the collagen IV molecules happens via LAMMPS
`fix bond/react` and a large set of template files. Each allowed reaction is
specified by three template files — a `pre` (structure of the bonding pattern
before the reaction), a `post` (after) and a `map` (instructions for
translating molecule identifiers and properties from `pre` to `post`).
Because of the many combinations of molecule types and valencies, there are
hundreds of template patterns. See
[`fix bond/react` in the LAMMPS docs](https://docs.lammps.org/fix_bond_react.html)
for the file format.

The two template families currently in the repo are:

| Folder | Reactions covered |
|---|---|
| `model_scripts/Stretch/templates/NC1_6N/` | NC1-valency change reactions (groups 01–15) |
| `model_scripts/Stretch/templates/OpenTetramerToRhomboidTetramer/` | Tetramer rearrangements (groups 01–41) |

> TODO: add a one-page index of which numerical group corresponds to which
> valency change.

---

## Data processing

`model_scripts/GaussianSpread/` converts LAMMPS trajectory frames into images
that can be analysed by the AFT (Alignment by Fourier Transform) tool
(see `readme.md` in that folder). This is the data-processing step used for
the analyses in Figures 3 and 4. The pipeline is:

1. `PointSpreadSingle.py` (or `_7s.py`) — render one trajectory frame as a
   point-spread image. `MetaScript.py` wraps this for many timepoints, and
   `GaussianSpreadImages.sh` runs that wrapper across many simulations.
2. `ImageCreation.py` — produces the final images / plots from the spread
   data.

`model_scripts/DataProcessing/` contains the other analyses:

- `MoleculeAngle.py` (`MoleculeAlignment.sh`) — molecular alignment
  (nematic order parameter) of the 7S beads.
- `Histogram_Types.py` (`HistogrammingTypes.sh`) — count molecule binding
  types over time.

These produce the `MolAlign*.txt` and related files that the plotting
notebooks consume.

---

## Reproducing the figures

| Figure | Folder | Notebook(s) | What it shows |
|--------|--------|-------------|---------------|
| Fig. 3 b–d | `Figures/Snapshots/Fig3GaussianSnapshots/` | `Simulation_Remodel_All.ovito` | Network self-assembly + stretch snapshots |
| Fig. 3 e–i | `Figures/FigAlignment/` | `Mechanics.ipynb`, `Alignment.ipynb` | Stress relaxation and nematic order under stretch, with vs without bond rearrangement |
| Fig. S3 (ablations) | `Figures/FigAblations/` | `PlottingAblations.ipynb` | Simulated annular ablations used to calibrate τ |
| Fig. 7 | `Figures/Fig7sMutants/` | `Mechanics7S.ipynb`, `Alignment.ipynb` | Networks with a fraction of low-valency 7S beads |
| Fig. 7 snapshots | `Figures/Snapshots/Fig7sMutantGaussian/` | OVITO scene | Mutant network snapshots before/after stretch |

`MergingFiles.ipynb` in `FigAlignment/` is a maintenance notebook — it
concatenates segments from successive `*RESTART` runs into single
`merged_seed*.txt` files. It does not produce a paper figure.

> TODO: confirm the figure numbering against the final paper version and add
> any other figures that are produced from this code.

---

## Data conventions

- **Seeds.** Replicate runs are labelled `seed1`, `seed2`, … in filenames.
- **Run conditions.** Folder names encode the simulation condition, e.g.
  - `NoStretch/` — equilibrated network, no stretch applied.
  - `NoRemodel_Xstretch100/` — uniaxial stretch of 1.0 with bond
    rearrangement disabled.
  - `Remodel_Xstretch100/` — same stretch, with bond rearrangement enabled.
  - `Equilibration_Fraction0.56/` — pre-stretch equilibration of a 7S-mutant
    network where 56% of 7S beads are wild-type valency.
- **File types.**
  - `thermo*.dat` — LAMMPS `thermo` output (energies, pressure tensor, box).
  - `MolAlignrun_seed*.txt` — molecular alignment (nematic order parameter)
    of the 7S beads as a function of time, produced by
    `model_scripts/DataProcessing/MoleculeAngle.py`.
  - `MolAlignrunRS_seed*.txt`, `MolAlignrunRSRS_seed*.txt`, … — the same
    quantity computed from successive RESTART segments, before they are
    merged by `MergingFiles.ipynb`.
  - `median_order_parameterWtNames_*.csv` — AFT-style order parameter
    measured from rendered snapshots (Gaussian-spread images), used as the
    simulation analogue of the experimental microscopy AFT.
  - `TimeSeries_Stretch*.csv`, `LastFrame_StretchMultiple.csv` — annular
    ablation observables.
- **Provenance.** `Figures/DataFrom.txt` records the original cluster paths
  the data in this repo was copied from; `Figures/Snapshots/FilesFrom.md`
  does the same for snapshot files.

> TODO: short table describing each column of `MolAlignrun_seed*.txt` and
> `thermo_seed*.dat` so users don't have to reverse-engineer the columns.

---

## Citing this work

If you use this code, please cite the paper:

```bibtex
@article{Khalilgharibi2026BM,
  title   = {Basement membrane sets the timescale of tissue mechanical memory},
  author  = {Khalilgharibi, Nargess and Meadowcroft, Billie and \v{S}ari\'c,
             An\d{e}la and Mao, Yanlan},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.64898/2026.04.23.720416}
}
```

---

## Contact

- Code questions: Billie Meadowcroft — billiemead@hotmail.co.uk
- Corresponding author: Yanlan Mao — y.mao@ucl.ac.uk

> TODO: add a licence file (MIT / BSD-3 / CC-BY-NC etc.) and reference it here.

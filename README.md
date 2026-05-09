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

## Table of contents

1. [Overview of the model](#overview-of-the-model)
2. [Repository layout](#repository-layout)
3. [Requirements](#requirements)
4. [Quick start](#quick-start)
5. [Reproducing the figures](#reproducing-the-figures)
6. [Data conventions](#data-conventions)
7. [Citing this work](#citing-this-work)
8. [Contact](#contact)

---

## Overview of the model

Each collagen IV protomer is represented as a two-bead rod: one bead for the
NC1 domain and one for the 7S domain. NC1 beads bind pairwise; 7S beads can
bind up to three other 7S beads, so protomers self-assemble into a percolating
network resembling the BM. Key simulation stages are:

1. **Self-assembly + equilibration** of the network from random initial conditions.
2. **Uniaxial stretch** of the equilibrated network at a chosen strain.
3. (Optional) **Bond rearrangement** — stochastic bond breaking and reforming
   without adding new protomers — to model network remodelling.
4. **Annular ablation** in silico, used to calibrate simulation timescales
   against the experimental ablation recoil.
5. **Mutant networks**: a fraction of 7S beads has reduced valency, mimicking
   the 7S mutation.

Outputs from the simulations are time series of stress, molecular alignment
(`MolAlign*`), order parameter snapshots (AFT analogue) and LAMMPS
trajectories.

> TODO: short reference to the specific equations / potentials, plus the
> previous model paper this builds on (ref. 37 in the manuscript).

---

## Repository layout

```
GithubCode/
├── README.md
├── lammps_src/          # LAMMPS input scripts + any custom source
│   └── TODO
├── model_scripts/       # Pre/post-processing helpers (initial configs, mutants, ablations)
│   └── TODO
└── Figures/             # One folder per figure in the paper
    ├── DataFrom.txt              # Provenance: where raw data on the cluster was generated
    ├── FigAlignment/             # Fig. 3 — alignment & stress, WT network
    │   ├── Alignment.ipynb
    │   ├── Mechanics.ipynb
    │   ├── MergingFiles.ipynb
    │   ├── Data/
    │   │   ├── Equilibration/
    │   │   ├── NoStretch/
    │   │   ├── NoRemodel_Xstretch100/
    │   │   └── Remodel_Xstretch100/
    │   └── Plots/
    ├── FigAblations/             # Fig. S3 / calibration — simulated annular cuts
    │   ├── PlottingAblations.ipynb
    │   └── Data_all/
    ├── Fig7sMutants/             # Fig. 7 — reduced 7S valency networks
    │   ├── Alignment.ipynb
    │   ├── Mechanics7S.ipynb
    │   └── Data/
    │       ├── Equilibration_Fraction0.56/
    │       ├── Equilibration_Fraction0.75/
    │       ├── XStretch100_Fraction0.56/
    │       └── XStretch100_Fraction0.75/
    └── Snapshots/                # Rendered network snapshots (OVITO scenes + tif/png)
        ├── FilesFrom.md
        ├── Simulation_Remodel_All.ovito
        ├── Fig3GaussianSnapshots/
        └── Fig7sMutantGaussian/
```

> TODO: once `lammps_src/` and `model_scripts/` are populated, fill in the
> filenames of the main run script, the ablation script, and the
> 7S-mutant generator.

---

## Requirements

**Simulations**

- LAMMPS (TODO: state version + required packages, e.g. MOLECULE, RIGID,
  USER-MISC, plus any custom fix used for bond breaking/forming)
- A C++ compiler if rebuilding LAMMPS with custom src files in `lammps_src/`
- MPI (optional) for parallel runs

**Analysis & plotting**

- Python ≥ 3.9
- `numpy`, `scipy`, `pandas`, `matplotlib`, `jupyter`
- AFT (Alignment by Fourier Transform) workflow — TODO: link / cite
- OVITO (for opening the `.ovito` snapshot scenes)

> TODO: add an `environment.yml` or `requirements.txt`.

---

## Quick start

```bash
# 1. clone
git clone <repo-url>
cd GithubCode

# 2. (TODO) build / install LAMMPS with the required packages
#    -- instructions in lammps_src/README

# 3. run a single equilibration + stretch trajectory
#    (TODO: example command, e.g.)
#    lmp -in lammps_src/run_stretch.in -var seed 1 -var xstretch 1.0

# 4. analyse and plot
jupyter lab Figures/FigAlignment/Alignment.ipynb
```

Each notebook expects to be run from the repo with the `Data/` folders sitting
next to it. The data paths inside the notebooks are currently hard-coded as
`/Users/billiemeadowcroft/Dropbox/Collagen/NargessPlotting/GithubCode/...`
— **TODO: replace these with relative paths so the notebooks are portable.**

---

## Reproducing the figures

| Figure | Folder | Notebook(s) | What it shows |
|--------|--------|-------------|---------------|
| Fig. 3 b–d | `Figures/Snapshots/Fig3GaussianSnapshots/` | `Simulation_Remodel_All.ovito` | Network self-assembly + stretch snapshots |
| Fig. 3 e–i | `Figures/FigAlignment/` | `Mechanics.ipynb`, `Alignment.ipynb` | Stress relaxation and nematic order under stretch, with vs without bond rearrangement |
| Fig. S3 (ablations) | `Figures/FigAblations/` | `PlottingAblations.ipynb` | Simulated annular ablations used to calibrate τ |
| Fig. 7 | `Figures/Fig7sMutants/` | `Mechanics7S.ipynb`, `Alignment.ipynb` | Networks with a fraction of low-valency 7S beads |
| Fig. 7 snapshots | `Figures/Snapshots/Fig7sMutantGaussian/` | OVITO scene | Mutant network snapshots before/after stretch |

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
    of the 7S beads as a function of time.
  - `median_order_parameterWtNames_*.csv` — AFT-style order parameter
    measured from rendered snapshots, used as the simulation analogue of the
    experimental microscopy AFT.
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

# Suggested improvements to the GithubCode repo

Goal: turn the repo into a clean, durable record of the model and simulations
behind the paper, so that a future reader (collaborator, reviewer, you in two
years) can (a) understand the model, (b) reproduce a figure, and (c) re-run a
simulation, without having to ask you anything.

The notes below assume **no code changes** — they are about structure,
documentation and hygiene. They are roughly ordered by impact.

---

## 1. High-impact "record of the paper" gaps

These are the things a reviewer or future user will notice first.

### 1.1 Add a LICENSE
There is no licence file. Without one, technically nobody is allowed to reuse
the code. Pick one (MIT or BSD-3 are conventional for academic simulation
code; CC-BY for the data/figures) and reference it in the README. The README
already has a TODO for this.

### 1.2 Add a CITATION.cff
GitHub renders this as a "Cite this repository" button. It takes about ten
lines and gives the repo a stable, machine-readable citation that's separate
from the paper citation. Useful when the bioRxiv DOI is superseded by the
journal DOI.

### 1.3 Pin software versions
Currently versioning is partial:
- LAMMPS version is named in `MyLooseInstructions.md` (`lammps-15Jun2023`)
  but **not** in the README. Move it into the README's *Requirements*
  section.
- No `environment.yml` / `requirements.txt`. Two notebooks already use
  `pandas`, `numpy`, `matplotlib`; `ImageCreation.py` also needs `seaborn`.
  Even an unpinned list is better than nothing.
- `model_scripts/Stretch/build_Stretch_SimpleAnnotated.py` documents the
  `bond/react` dependence; mention which LAMMPS package(s) need to be
  enabled at build time (`MOLECULE`, `RIGID`, `REACTION`).

### 1.4 Resolve the README TODOs
The README is well-structured but has eight `TODO` markers. The most
load-bearing ones:
- `lammps_src/` and `model_scripts/` "TODO" placeholders in the *Repository
  layout* tree. These now contain real files — update the tree to show what
  is actually there.
- "TODO: example command" in *Quick start*. Even a single one-line
  invocation that you've actually run is enormously useful. Without it the
  Quick start section is essentially empty.
- "TODO: short table describing each column of `MolAlignrun_seed*.txt` and
  `thermo_seed*.dat`". This is the single highest-leverage doc to add.
  Without it, anyone re-using the data has to reverse-engineer the columns
  from the notebooks.

### 1.5 Merge `MyLooseInstructions.md` into the rest of the docs
Right now there are two parallel sets of instructions:
- `README.md` — polished, paper-facing.
- `MyLooseInstructions.md` — practical, has the actually-tested commands and
  cluster paths.

The "loose" file is full of valuable knowledge (LAMMPS version, exact
`mpirun` invocation, equilibration step counts, why
`build_Stretch_SimpleAnnotated.py` exists, why type 9 is needed for 7s
mutants, the "TO do (myself)" list). Promote the user-facing parts into the
README (under *Quick start*, *Reproducing the figures*, and a new *Model
details* section), and either delete `MyLooseInstructions.md` or rename it
to `DEVELOPMENT_NOTES.md` and explicitly mark it as personal scratch.

The "TO do (myself)" list at the bottom should also be tracked as GitHub
issues, not in a markdown file.

---

## 2. Make the repo runnable end-to-end

The repo currently *describes* the workflow but doesn't quite let someone
follow it. A few small additions close the gap.

### 2.1 Per-stage READMEs
The pattern `model_scripts/Initial_config/Readme.md` and
`model_scripts/GaussianSpread/readme.md` is good. Add one per stage:
- `model_scripts/Stretch/README.md` — what `Stretch.sh`, `Stretch_NoRemodel.sh`
  and the `*RESTART.sh` variants do, in what order, and which `build_*.py` is
  called by which shell wrapper.
- `model_scripts/7s_Perturbations/README.md` — same, plus the special role of
  bead type 9.
- `model_scripts/Ablations/README.md` — explain
  `SimulationsStartFromRestartFile/run_data4_T1.465000.restart` (where it came
  from, which equilibration produced it, what its box size is).
- `model_scripts/DataProcessing/README.md` — what each `*.sh` / `*.py`
  produces, and which figure(s) consume that output.
- `lammps_src/README.md` — that LAMMPS version `15Jun2023` must be patched
  by replacing `src/REACTION/fix_bond_react.{cpp,h}` with the files here,
  what the modification does, and why it was needed.

These are short — half a page each — but together they replace the
"figure out which script feeds which" investigation that a new user
otherwise has to do.

### 2.2 Tighten the existing per-stage READMEs
- `Initial_config/Readme.md` is a notes file (`Bond lenggth`, three different
  scenarios listed without context). Rewrite it as: "to generate config X,
  run `python create_dumbbell_lattice_slabNC1_type9.py` and answer the prompts
  with [these values]; the result is `dumbbells_*.lammpsdata` which is what
  `Equilibrate.sh` reads."
- `GaussianSpread/readme.md` is full of personal cluster paths. Replace
  `/nfs/scistore26/saricgrp/bmeadowc/Scratch/...` with a placeholder
  (e.g. `<RUN_DIR>`) and explain the layout the script expects.

### 2.3 Add an end-to-end "minimal example"
A single script (or section in the top-level README) that runs *one* short
seed of *one* condition and produces *one* figure panel, in well under an
hour, is the gold standard for paper-companion repos. It validates the whole
pipeline at once and gives users something to diff against. Even just
documenting "this is the smallest meaningful run" without scripting it would
help.

### 2.4 Snapshot a small example dataset for the notebooks
The notebooks already have full `Data/` folders alongside them — good. Make
that explicit in each notebook (a markdown cell at the top: "this notebook
reads from `./Data/`, which is bundled with the repo. To regenerate the
data, see `model_scripts/Stretch/README.md`").

---

## 3. Notebook hygiene

Every plotting notebook starts with absolute paths like:

```python
datadir0  = '/Users/billiemeadowcroft/Dropbox/Collagen/NargessPlotting/GithubCode/Figures/FigAlignment/Data/'
plotsdir  = '//Users/billiemeadowcroft/Dropbox/Collagen/NargessPlotting/GithubCode/Figures/FigAlignment/Plots/'
```

(note the doubled `//` — a real bug waiting to happen).

This is the single biggest barrier to anyone else running the notebooks. The
README already flags it as a TODO. Consistent recommendation:

- Replace with `Path(__file__).parent / "Data"` or, for notebooks,
  `Path.cwd() / "Data"` after a comment saying "run from this notebook's
  directory".
- Move the colour palette (`clrs = [...]`) and constants like
  `NumMolecules = 15680`, `TEq = 1085000`, `Vol = 164*164*12` to a tiny
  shared `figure_style.py` / `constants.py`. They are repeated across
  `Mechanics.ipynb`, `Alignment.ipynb`, `Mechanics7S.ipynb`, etc.
- "Clear all outputs" before committing notebooks, or use
  [`nbstripout`](https://github.com/kynan/nbstripout). Stored outputs make
  diffs unreadable and bloat the repo.
- Each notebook should open with a one-paragraph markdown cell stating which
  figure(s) it produces, what data it expects in `Data/`, and what files it
  writes into `Plots/`.

Also: `Figures/FigAlignment/MergingFiles.ipynb` is clearly a maintenance
notebook (concatenates restart segments). Move it to a `tools/` subfolder
or rename it `_MergingFiles.ipynb` so it's obvious it's not figure-producing.

---

## 4. Repository hygiene

Easy, mechanical fixes:

### 4.1 Add a `.gitignore`
Currently the repo tracks
- `.DS_Store` files (macOS Finder metadata) — at least four of them.
- `model_scripts/Stretch/templates/OpenTetramerToRhomboidTetramer/some_files_changes.ipynb` — name suggests scratch.
- `Figures/FigAlignment/Data/Remodel_Xstretch100/MergingFiles/` contains
  `MolAlignrunRSRSRS_seed*.txt`, `merged_seed*.txt`, `MolAlignrun_seed4o.txt` —
  intermediate products that the notebook produces. Either commit only the
  final merged outputs, or commit only the raw inputs and document the
  merging step.
- `Figures/FigAlignment/Data/NoStretch/oldthermo/` and
  `MolAlignrun_seed*_o.txt` (the `_o` / `oldthermo` suffix suggests
  superseded data). Either delete or label clearly.

A minimal `.gitignore`:
```
.DS_Store
__pycache__/
*.pyc
.ipynb_checkpoints/
```
plus optionally `*.restart` and `*.lammpsdata` if you don't want large binary
inputs in git (see 4.2).

### 4.2 Big binary files — consider Git LFS
The repo currently contains `.tif` snapshots, `.ovito` scenes, a `.restart`
file, a `.lammpsdata` initial config and `.dat` thermo files. None of these
are huge individually but together they bulk up clones, and they don't diff
usefully. Two reasonable policies:
- "Source-only" repo: move binaries to a Zenodo / OSF deposit and link from
  the README. Add a small `Figures/Snapshots/README.md` explaining how to
  regenerate or download them.
- "Self-contained" repo: keep them but use [Git LFS](https://git-lfs.com/)
  for `*.tif`, `*.ovito`, `*.restart`, `*.lammpsdata`, `*.lammpstrj`.

Either is fine; the current setup (raw binaries in main git) becomes
painful as the repo grows.

### 4.3 Naming consistency
- Top-level folder is `GithubCode/`. When the repo lives at
  `github.com/<org>/<repo>` this becomes redundant. Consider renaming the
  GitHub repo itself something like `bm-mechanical-memory-sims` and let the
  inner folder be the repo root.
- Mixed casing: `Readme.md` (Initial_config), `readme.md` (GaussianSpread),
  `README.md` (root). Pick one — `README.md` is conventional and is the
  only one GitHub auto-renders inside subfolders.
- `model_scripts/` is `snake_case` in spirit but `Initial_config`,
  `GaussianSpread`, `DataProcessing`, `Ablations` are `PascalCase`. Same
  for `Figures/Fig*`. Not worth a churn-only PR, but if you ever rename
  things, pick one convention.

### 4.4 The `templates/NC1_6N/` directory (~hundreds of files)
This is unavoidable — `bond/react` requires one `pre/post/map` triple per
allowed reaction. But:
- Add a `templates/README.md` that explains the naming convention
  (`NN_pre.txt`, `NN_post.txt`, `NN_map.txt`, `*_frag.txt`,
  `NN_map_break.txt`, `NN_map_creation.txt`) and what each numerical group
  (01–15) represents (which valency change).
- A small index table: "reaction 03 = NC1 with 2 bonds + free NC1 → NC1
  with 3 bonds + bonded NC1" etc. If templates were generated by a script,
  commit that script and link to the LAMMPS doc you already mention
  (`https://docs.lammps.org/fix_bond_react.html`).

---

## 5. Provenance of data

This is the part that decides whether the repo is a useful *record* of the
paper or just a convenient toolbox.

### 5.1 Strengthen the cluster-path provenance
`Figures/DataFrom.txt` and `Figures/Snapshots/FilesFrom.md` are good ideas
— but they're inconsistent (one is `.txt`, one is `.md`), partial (no
provenance for `Fig7sMutants/Data` or `FigAblations/Data_all`), and they
list ISTA cluster paths that won't survive your account being archived.

Suggestions:
- One provenance file per `Data/` folder (`Data/PROVENANCE.md`) with: the
  cluster path the data came from; the git commit of `model_scripts/`
  used to produce it; the parameter set (strain, density, seed range);
  and the date.
- Or, equivalently, embed the parameters directly in the *folder name*
  (you're already doing this with `XStretch100_Fraction0.56` etc., which
  is great) plus a single `PROVENANCE.md` per figure folder summarising
  all conditions.

### 5.2 Document what each data file column means
The single most important table to add. Suggested format, in
`Figures/README.md`:

| Filename pattern | Producer | Columns | Units |
|---|---|---|---|
| `thermo_seed*.dat` | LAMMPS thermo output via `build_Stretch.py` | step, temp, pe, ke, stressX, stressY, stressZ, ... | LJ |
| `MolAlignrun_seed*.txt` | `model_scripts/DataProcessing/MoleculeAngle.py` | per-molecule cos(θ) values, one frame per block | dimensionless |
| `median_order_parameterWtNames_*.csv` | AFT pipeline on Gaussian-spread images | timestep, median_OP | dimensionless |
| `TimeSeries_Stretch*.csv` | `PlottingAblations.ipynb` upstream | … | … |
| `LastFrame_StretchMultiple.csv` | … | … | … |

Even partial completion of this table is a huge win.

### 5.3 Capture parameters in machine-readable form
At the moment, simulation parameters live inside shell scripts as bash
variables (`Time_Str=2e5`, `BP=0.01`, …). Consider a small
`parameters.yaml` per simulation condition (committed alongside the data),
so a notebook can `yaml.safe_load()` them rather than hard-coding e.g.
`TEq = 1085000` in three notebooks. This is the kind of subtle change that
pays off the next time someone wants to run a slightly-different parameter
sweep.

---

## 6. Smaller polish items

- README `Overview of the model` reads well; the one missing piece is a
  pointer to the previous paper this builds on (currently a TODO marked as
  "ref. 37 in the manuscript"). A direct DOI in the README means the
  reader doesn't have to dig.
- Snapshots: `Figures/Snapshots/Fig7sMutantGaussian/` contains files named
  `Time2727000_PostStretch_seed2__Frac0.56Mutant.tif` (note the double
  underscore). Worth normalising filenames if you do another pass.
- `Figures/FigAblations/PlottingAblations.ipynb` — confirm it actually
  produces Fig. S3 panels and add the figure number(s) in a top markdown
  cell.
- `build_Stretch_SimpleAnnotated.py` is genuinely excellent — that
  pattern of "annotated, non-runnable companion to the real script"
  could be extended to `build_StretchAblate.py` and `build_Equilibrate.py`,
  which are similarly long and template-heavy.
- The fact that `build_Stretch_SimpleAnnotated.py` exists is itself the
  best feature of the repo for "record of the paper" purposes; mention
  it prominently in the top-level README under *Overview of the model*.

---

## 7. Suggested README rewrite skeleton

After the above, the top-level README would look like:

1. Title + paper citation (already done).
2. Overview of the model — short, with link to the predecessor paper.
3. Repository layout — updated tree, with one-line descriptions, no TODOs.
4. Requirements — LAMMPS version + packages, Python `requirements.txt`,
   AFT pipeline link, OVITO.
5. Quick start — *one* worked example, end-to-end, with copy-pasteable
   commands.
6. Reproducing each figure — table you already have, but with a column for
   "expected runtime" so users know if it's a 1-minute notebook or a
   500-CPU-hour simulation.
7. Data conventions — promoted version of the column-meaning table.
8. Citing this work + LICENSE + CITATION.cff.
9. Contact.

`MyLooseInstructions.md` disappears (its content is now distributed where
it belongs); per-stage READMEs replace the long shell-command-by-shell-command
narrative.

---

## Priority order if you only do five things

1. Add a LICENSE.
2. Document the columns of `thermo_seed*.dat` and `MolAlignrun_seed*.txt`.
3. Replace absolute paths in notebooks with relative paths (the doubled
   `//` is a latent bug too).
4. Add `.gitignore` + clear out `.DS_Store`, `oldthermo/`, `*_o.txt`,
   `MergingFiles/` intermediates.
5. Promote the useful parts of `MyLooseInstructions.md` into the README
   and per-stage READMEs; delete the rest.

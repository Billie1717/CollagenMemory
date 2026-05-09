"""
build_Stretch_SimpleAnnotated.py
================================

Simplified, *annotated* version of `build_Stretch.py`. This script is for
**reading**, not for running -- the working version is `build_Stretch.py`,
which is what the paper actually uses.

What `build_Stretch.py` does:
    1. Reads simulation parameters from the command line.
    2. Writes a LAMMPS input script (`collagen.in`) that:
         - reads an equilibrated network from a restart file,
         - sets pair / bond / angle styles and force-field constants,
         - declares hundreds of `molecule` templates and `react` rules
           (one pair per allowed change of NC1 valency), and uses the
           `bond/react` package to form / break bonds stochastically,
         - applies a uniaxial stretch via `fix deform`,
         - relaxes the network for `trelax` steps with bond rearrangement
           on (the "Remodel" condition).
    3. Writes a SLURM submission script (`runscript.sh`) that launches
       LAMMPS on a single core.

In this simplified version we keep one example of each kind of LAMMPS
command and replace the long lists of templates / reactions with a single
representative line plus a comment explaining the pattern. Anywhere a
`...` appears in the LAMMPS strings means "in the real script there are
many more lines like the one above".

The script is therefore NOT runnable -- it's a guided tour of the input
file. Read it top-to-bottom alongside `build_Stretch.py` to understand
what each block of the real script does.
"""

from __future__ import division, print_function
import sys


def main():
    # ------------------------------------------------------------------
    # 1. Command-line parameters
    # ------------------------------------------------------------------
    # The launcher script `Stretch.sh` passes these in as positional
    # arguments. Order matters and is the same as in `build_Stretch.py`.
    workingdir = str(sys.argv[1])    # where collagen.in / runscript.sh land
    tStr       = str(sys.argv[2])    # number of steps over which to apply stretch
    tRel       = str(sys.argv[3])    # number of steps to relax after stretch
    tInitRel   = str(sys.argv[4])    # short relaxation between stretch and remodelling
    Frames     = str(sys.argv[5])    # number of trajectory dumps requested
    Nevery     = str(sys.argv[6])    # how often (in steps) bond/react is attempted
    tstep      = str(sys.argv[7])    # MD integration timestep (LJ units)
    MakeDist   = str(sys.argv[8])    # max distance at which a bond can form (RmaxF)
    BreakDist  = str(sys.argv[9])    # min distance at which a bond can break (RminFB)
    MakeProb   = str(sys.argv[10])   # per-attempt probability of forming a bond
    BreakProb  = str(sys.argv[11])   # per-attempt probability of breaking a bond
    Xstretch   = str(sys.argv[12])   # uniaxial strain as a fraction of box length
    seed       = str(sys.argv[13])   # RNG seed (used for filenames / metadata)

    # ------------------------------------------------------------------
    # 2. Build collagen.in
    # ------------------------------------------------------------------
    # The real script writes ~1100 lines into collagen.in by concatenating
    # many `fo.write(...)` calls. We do the same here, but each block is
    # cut down to the essential commands and heavily commented.
    print("Writing collagen.in (annotated demo)\n")
    fo = open(workingdir + "/collagen.in", "w")

    # ---- 2a. Run-length variables ------------------------------------
    # Every `variable name equal value` line just defines a LAMMPS
    # variable, so the input file is parameterised rather than hard-coded.
    fo.write(f"""\
#####################################################
# RUN LENGTHS (in MD steps)
variable tinitRel  equal {tInitRel}   # short relax between stretch and remodelling
variable tstretch  equal {tStr}       # steps over which the stretch is ramped
variable trelax    equal {tRel}       # main relaxation run with bond/react ON
variable Frames    equal {Frames}     # how many trajectory frames we want

# Derived dump periods (so output is independent of run length)
variable ttotal     equal $(v_trelax+v_tstretch)
variable thermodump equal $(floor(v_ttotal/5e4))   # thermo line every N steps
variable trajdump   equal $(floor(v_ttotal/v_Frames))
variable trestart   equal $(floor(v_ttotal/20))

#####################################################
# THERMOSTAT / INTEGRATOR PARAMETERS
variable mytemp equal 1.0           # Langevin temperature (LJ)
variable damp   equal 0.1           # Langevin damping coefficient
variable tstep  equal {tstep}       # MD timestep
variable lj_cutoff equal 1.122462   # WCA cutoff = 2^(1/6)

# ---- bond/react parameters (passed in from the launcher) ----
variable RmaxF   equal {MakeDist}    # form a bond if pair within this distance
variable probF   equal {MakeProb}    # per-attempt formation probability
variable RminFB  equal {BreakDist}   # break a bond beyond this distance
variable probFB  equal {BreakProb}   # per-attempt breaking probability
variable RminF   equal 0.0           # lower form cutoff (always 0)
variable RmaxFB  equal 30            # large outer cutoff for break neighbour search
variable NeverySlow equal {Nevery}   # bond/react attempt frequency

# RNG seeds for the 14 different reaction families. In the real script these
# are 14 separate variables (lseed1 .. lseed14) so each family can have an
# independent stream; here we just show the pattern.
variable lseed_base equal 13579
variable lseed1     equal ${{lseed_base}}+1   # NC1 creation
# variable lseed2 ... lseed14   <-- one per reaction family in build_Stretch.py
""")

    # ---- 2b. Force-field constants -----------------------------------
    fo.write("""
#####################################################
# FORCE FIELD
units        lj
boundary     p p f                 # periodic in x,y; fixed walls in z
atom_style   molecular
bond_style   harmonic
angle_style  harmonic
pair_style   hybrid/overlay lj/cut 1.0 zero 1.0
pair_coeff   * * lj/cut 0.0 0.0 ${lj_cutoff}   # WCA only -- purely repulsive
pair_coeff   * * zero  6.0                     # placeholder for neighbour list

# Spring constants (k) and rest lengths (r0)
variable k_nc1  equal 6.0
variable r0_nc1 equal 0.5
variable k_7s   equal 6.0
variable r0_7s  equal 0.5

bond_coeff  1  100.0 3.0           # rigid backbone bond inside a protomer
bond_coeff  2  ${k_nc1} ${r0_nc1}  # NC1-NC1 bond (between protomer ends)
bond_coeff  3  ${k_7s}  ${r0_7s}   # 7S-7S bond (multi-valent end)

angle_coeff 1  1.0 180.0
angle_coeff 2  1.0 155.0
angle_coeff 3  1.0  60.0

variable max_bonds_7s equal 3      # max valency at the 7S bead
""")

    # ---- 2c. Read the equilibrated configuration ---------------------
    fo.write("""
#####################################################
# INITIAL STATE
# We start from a restart file produced by `build_Equilibrate.py`, which
# already contains the self-assembled network.
variable fname     string data
variable simname   string run_T${mytemp}
read_restart       ${fname}

variable box_size       equal xhi-xlo
variable x_stretch_Frac equal """ + Xstretch + """
variable x_stretch      equal $(v_x_stretch_Frac*v_box_size/2)

neighbor      0.4 bin
neigh_modify  every 10 delay 0 check yes
timestep      ${tstep}
""")

    # ---- 2d. Molecule templates --------------------------------------
    # In `build_Stretch.py` this block is enormous: every allowed local
    # bonding topology (NC1 with 2,3,...,8 neighbours, plus loops, dimers,
    # line/triangle trimers, open/rhomboid/square tetramers, etc.) has a
    # `_pre.txt` template (state before the reaction) and a `_post.txt`
    # template (state after). We show one example below.
    fo.write("""
#####################################################
# MOLECULE TEMPLATES  (~700 lines in the full script; one example shown)
# Each reaction needs a "pre" template (matching graph) and a "post"
# template (graph after the bond change). Filenames live under
# ../../templates/<FamilyName>/<ID>_(pre|post|map).txt
molecule  Nc1_2N_01_pre  ../../templates/NC1_2N/01_pre.txt
molecule  Nc1_2N_01_post ../../templates/NC1_2N/01_post.txt
# ... Nc1_2N_02 .. Nc1_8N_01 (NC1 valency 2..8)
# ... Nc1_loop_*           (closed loops)
# ... MonoToDimer_*        (n-mer growth: monomer -> dimer)
# ... DimerToLineTrimer_*, LineTrimerToTriangleTrimer_*
# ... TriangleTrimerToOpenTetramer_*, OpenTetramerToRhomboidTetramer_*
# ... RhomboidTetramerToSquareTetramer_*  (and all reverse moves)
""")

    # ---- 2e. Computes / thermo / dumps -------------------------------
    fo.write("""
#####################################################
# OBSERVABLES
compute  perAtomStress all stress/atom NULL
compute  totalStressX  all reduce sum c_perAtomStress[1]
compute  totalStressY  all reduce sum c_perAtomStress[2]
compute  totalStressZ  all reduce sum c_perAtomStress[3]
compute  forceBond     all bond/local force
compute  forceBondAll  all reduce ave c_forceBond

thermo_style  custom step etotal ke pe temp press &
              c_totalStressX c_totalStressY c_totalStressZ
thermo        ${thermodump}

# Trajectory + bond list dumps (one frame every trajdump steps)
dump  2 all custom ${trajdump} dumplin/dump.${simname}.*.lammpstrj &
        type id xu yu zu
dump  3 all local  ${trajdump} dumplin_bonds/bonds.${simname}.* &
        index c_1[1] c_1[2] c_1[3]
restart      ${trestart} restart/${simname}.*.restart
""")

    # ---- 2f. Thermostat + initial velocities -------------------------
    fo.write("""
#####################################################
# NVE + LANGEVIN
fix  fix_wall_lo  all wall/lj126 zlo EDGE 1.0 1.0 ${lj_cutoff}
fix  fix_wall_hi  all wall/lj126 zhi EDGE 1.0 1.0 ${lj_cutoff}
fix  fix_lengevin all langevin  ${mytemp} ${mytemp} ${damp} ${lseed1}
fix  fNVE         all nve

velocity all create ${mytemp} 12345 dist gaussian mom yes rot yes
""")

    # ---- 2g. The actual experiment: stretch, settle, then remodel ----
    # This is the heart of the script -- everything above is setup.
    fo.write("""
#####################################################
# EXPERIMENT
write_restart restart.postmakebonds

# (1) Apply uniaxial stretch over `tstretch` steps. `fix deform ... remap x`
#     ramps the x-extent of the simulation box and remaps atom positions
#     accordingly, producing a quasi-static stretch.
fix   fDeform all deform 1 x delta -$(v_x_stretch) $(v_x_stretch) remap x
run   ${tstretch}
unfix fDeform
write_restart restart.poststretch

# (2) Short relaxation with NO bond rearrangement (network adjusts elastically).
run   ${tinitRel}

# (3) Long relaxation WITH bond rearrangement -- this is the remodelling
#     phase that build_StretchNoRemodel.py omits. The `bond/react` fix
#     attempts every reaction every NeverySlow steps and accepts it with
#     probability probF (forming) or probFB (breaking). In the real input
#     file this fix is huge: ~hundreds of `react ...` lines, one per
#     template pair. The two examples below show the structure:
fix freact_creation all bond/react stabilization no reset_mol_ids no &
    # forward reactions (form a bond) -- use *_pre -> *_post
    react Rc1_2N_01 all ${NeverySlow} ${RminF} ${RmaxF} &
          Nc1_2N_01_pre Nc1_2N_01_post &
          ../../templates/NC1_2N/01_map.txt prob v_probF ${lseed1} &
    # reverse reactions (break a bond) -- swap pre/post and use the
    # break probability + outer cutoffs
    react BRc1_2N_01 all ${NeverySlow} ${RminFB} ${RmaxFB} &
          Nc1_2N_01_post Nc1_2N_01_pre &
          ../../templates/NC1_2N/01_map.txt prob v_probFB ${lseed1}
    # ... ~250 more `react` lines covering every (Family_ID, direction)
    #     pair declared in the molecule-templates block above.

run ${trelax}
""")

    fo.close()

    # ------------------------------------------------------------------
    # 3. Build runscript.sh (SLURM submission for the cluster)
    # ------------------------------------------------------------------
    # The real script is identical in spirit -- just SBATCH headers plus
    # an `mpirun lmp_mpi -in collagen.in` line. We strip the headers down
    # to the ones that actually matter for understanding.
    print("Writing runscript.sh (annotated demo)\n")
    f = open(workingdir + "/runscript.sh", "w")
    f.write("""#!/bin/bash
#SBATCH -J Stretch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=100:00:00

module load gcc openmpi

# bond/react needs RNG seeds in two files; LAMMPS reads them via the
# `variable lseed file lseed.dat` line in collagen.in.
mkdir dumplin dumplin_bonds restart
printf $RANDOM > lseed.dat
printf $RANDOM > vseed.dat

# Run LAMMPS (must be the custom build that includes the modified
# REACTION package -- see lammps_src/ in the repo).
mpirun -np 1 lmp_mpi -in collagen.in
""")
    f.close()


if __name__ == "__main__":
    main()

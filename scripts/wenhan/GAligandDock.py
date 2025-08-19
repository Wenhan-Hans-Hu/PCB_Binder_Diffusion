## Wenhan: This script utilizes GALigandDock primarily for energy calculations. 
## To serve this purpose, I have selected runmode = 'VSH':
'''
Virtual screening-High accracy Allow receptor flexibility and more rigorous entropy calculation. Single run recommended for efficiency. Each takes 10~15 minutes.
<GALigandDock name="dock" runmode="VSH" scorefxn="dockscore" scorefxn_relax="relaxscore" nativepdb="holo.pdb"/>
'''
## The description above is for reference only, please DO NOT attempt to use it directly, as it applies to Rosetta, not PyRosetta.
## Requirements to run:
##   1) A ligand mol2 file (ligand.mol2)
##   2) A rosetta ligand.param file (LIG.params)
##   3) A PDB file containing the ligand in the protein pocket (holo.pdb) 
## Note: Make sure contents of your ligand of interest locate at the end of holo.pdb, which means HETATM lines for any cofactor(s) / metals / waters should locate above it.
## Note: No need to take care about conformational information. Its will treat your ligand as rotamer (which is inevitable) and score for energy.
## Note: Use diffusion environment or add pyrosetta in pipeline environment
## Written on 27 Mar 2025

## 04 Apr 2025 Update: Dock Rigid runmode is selected to use instead of `VSH`, as VSH require .cst file to apply necessary restraints:
'''
Self docking Self docking, when no change in receptor conformation is expected. 1~5 repeats recommended, each takes 3~10 minutes:
<GALigandDock name="dock" runmode="dockrigid" scorefxn="dockscore" scorefxn_relax="relaxscore" />
'''

###        Please customise this PBS generation script before use:
'''
### Setting up design directory and commands
os.chdir(WDIR)
Score_GAligandDock = f"{WDIR}/Score_GAligandDock_Directory"
os.makedirs(Score_GAligandDock, exist_ok=True)
os.chdir(Score_GAligandDock)

holo_pdb_DIR = f"{WDIR}/input_directory"
os.makedirs(Score_GAligandDock+"/logs", exist_ok=True)

commands_GAligandDock = []
cmds_filename_des = "commands_GAligandDock"
with open(cmds_filename_des, "w") as file:
    for pdb in glob.glob(f"{holo_pdb_DIR}/*.pdb"):
        commands_GAligandDock.append(f"{PYTHON['general']} {SCRIPT_DIR}/scripts/design/GAligandDock.py "
                             f"--pdb {pdb} --nstruct {NSTRUCT} "
                             f"--params {' '.join(params)} > logs/{os.path.basename(pdb).replace('.pdb', '.log')}\n") 
        file.write(commands_design[-1])

print("Example design command:")
print(commands_design[-1])


### Running design jobs with Slurm.
submit_script = "submit_GAligandDock.sh"
utils.create_slurm_submit_script(filename=submit_script, name="GAligandDock", mem="4g", gpu=True, group=50, 
                                 N_cores=1, time="3:00:00", email=EMAIL, array=len(commands_design),
                                 array_commandfile=cmds_filename_des)
'''
import os, argparse, time, sys
import pyrosetta as pyr
import pyrosetta.rosetta
import pandas as pd
from pyrosetta.rosetta.protocols.ligand_docking.ga_ligand_dock import GALigandDock
from pyrosetta.rosetta.core.scoring import ScoreType, n_score_types

# -------- Argument Parsing -------- #
parser = argparse.ArgumentParser()
parser.add_argument("--pdb", required=True, type=str, help="Input PDB file")
parser.add_argument("--params", nargs="+", type=str, help=".params file(s) for ligand/cofactors")
parser.add_argument("--nstruct", type=int, default=1, help="Number of independent docking runs")
parser.add_argument("--outdir", type=str, default="outputs", help="Directory for output PDBs")
parser.add_argument("--runmode", type=str, default="dockrigid", choices=["dockrigid", "dockflex", "VSH", "VSX"],
                    help="GALigandDock runmode")
args = parser.parse_args()

# -------- Setup -------- #
if not os.path.exists(args.outdir):
    os.makedirs(args.outdir, exist_ok=True)

extra_res_fa = ""
if args.params:
    extra_res_fa = "-extra_res_fa " + " ".join(args.params)

pyr.init(f"{extra_res_fa} "
         "-corrections:gen_potential "
         "-score::hb_don_strength hbdon_GENERIC_SC:1.45 "
         "-score::hb_acc_strength hbacc_GENERIC_SP2SC:1.19 "
         "-score::hb_acc_strength hbacc_GENERIC_SP3SC:1.19 "
         "-score::hb_acc_strength hbacc_GENERIC_RINGSC:1.19 "
         "-no_autogen_cart_improper")

pose = pyrosetta.pose_from_file(args.pdb)
Ligand_resno = pose.size()
if not pose.residue(Ligand_resno).is_ligand():
    raise ValueError("Last residue is not a ligand. Please check input PDB.")
ligand_chain_id = pose.pdb_info().chain(Ligand_resno)
print(f"Detected ligand chain: {ligand_chain_id}")

dockscore = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function("beta_genpot")
dockscore.set_weight(pyrosetta.rosetta.core.scoring.fa_rep, 0.2)
dockscore.set_weight(pyrosetta.rosetta.core.scoring.coordinate_constraint, 0.1)

dock = GALigandDock()
try:
    dock.setup_params_for_runmode(args.runmode)
    print(f"GALigandDock runmode '{args.runmode}' configured.")
except Exception as e:
    print(f"[!] Failed to configure runmode '{args.runmode}': {e}")
    sys.exit(1)

scorefile = os.path.join(args.outdir, "Scoring.csv")
all_scores = []

# -------- Docking Loop -------- #
for N in range(args.nstruct):
    print(f"\n=== Starting docking run {N + 1}/{args.nstruct} ===")
    t0 = time.time()

    working_pose = pose.clone()

    try:
        dock.apply(working_pose)
    except Exception as e:
        print(f"[!] GALigandDock failed on run {N}: {e}")
        continue

    output_pdb = os.path.join(args.outdir, f"{os.path.splitext(os.path.basename(args.pdb))[0]}_{N}.pdb")
    working_pose.dump_pdb(output_pdb)
    print(f"Saved docked pose to: {output_pdb}")

    # Initialize score_data with minimal scores
    total = dockscore(working_pose)
    score_data = {
        "run_id": N,
        "output_name": os.path.basename(output_pdb),
        "total_score": total,
        "score_per_res": total / working_pose.size()
    }

    # Attempt to extract energy breakdown
    try:
        emap = working_pose.energies().total_energies()
        for i in range(1, n_score_types):  # start from 1 to skip invalid 0
            stype = ScoreType(i)
            value = emap[stype]
            if abs(value) > 1e-6:
                score_data[str(stype)] = value
    except Exception as e:
        print(f"[!] Failed to extract energy terms on run {N}: {e}")

    all_scores.append(score_data)
    print(f" Run {N} finished in {round(time.time() - t0, 2)} seconds.")

# -------- Append to Score CSV -------- #
if all_scores:
    df = pd.DataFrame(all_scores)
    write_header = not os.path.exists(scorefile)
    df.to_csv(scorefile, mode='a', header=write_header, index=False)
    print(f"\nAppended scores to: {scorefile}")
else:
    print("\nNo successful runs. CSV not updated.")

'''
import os, argparse, time, sys
import pyrosetta as pyr
import pyrosetta.rosetta
import pandas as pd
from pyrosetta.rosetta.protocols.ligand_docking.ga_ligand_dock import GALigandDock
parser = argparse.ArgumentParser()
parser.add_argument("--pdb", required=True, type=str, help="Input PDB file")
parser.add_argument("--params", nargs="+", type=str, help=".params file(s) for ligand/cofactors")
parser.add_argument("--nstruct", type=int, default=1, help="Number of independent docking runs")
parser.add_argument("--outdir", type=str, default="outputs", help="Directory for output PDBs")
parser.add_argument("--runmode", type=str, default="dockrigid", choices=["dockrigid", "dockflex", "VSH", "VSX"],
                    help="GALigandDock runmode")
args = parser.parse_args()
if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

extra_res_fa = ""
if args.params:
    extra_res_fa = "-extra_res_fa " + " ".join(args.params)

pyr.init(f"{extra_res_fa} "
         "-corrections:gen_potential "
         "-score::hb_don_strength hbdon_GENERIC_SC:1.45 "
         "-score::hb_acc_strength hbacc_GENERIC_SP2SC:1.19 "
         "-score::hb_acc_strength hbacc_GENERIC_SP3SC:1.19 "
         "-score::hb_acc_strength hbacc_GENERIC_RINGSC:1.19 "
         "-no_autogen_cart_improper")

pose = pyrosetta.pose_from_file(args.pdb)
Ligand_resno = pose.size()
if not pose.residue(Ligand_resno).is_ligand():
    raise ValueError("Last residue is not a ligand. Please check input PDB.")
ligand_chain_id = pose.pdb_info().chain(Ligand_resno)
print(f"Detected ligand chain: {ligand_chain_id}")

dockscore = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function("beta_genpot")
dockscore.set_weight(pyrosetta.rosetta.core.scoring.fa_rep, 0.2)
dockscore.set_weight(pyrosetta.rosetta.core.scoring.coordinate_constraint, 0.1)

dock = GALigandDock()
try:
    dock.setup_params_for_runmode(args.runmode)
    print(f"GALigandDock runmode '{args.runmode}' configured.")
except Exception as e:
    print(f"[!] Failed to configure runmode '{args.runmode}': {e}")
    sys.exit(1)
scorefile = os.path.join(args.outdir, "Scoring.csv")
all_scores = []

for N in range(args.nstruct):
    print(f"\n=== Starting docking run {N + 1}/{args.nstruct} ===")
    t0 = time.time()

    working_pose = pose.clone()

    try:
        dock.apply(working_pose)
    except Exception as e:
        print(f"GALigandDock failed on run {N}: {e}")
        continue

    output_pdb = os.path.join(args.outdir, f"{os.path.splitext(os.path.basename(args.pdb))[0]}_{N}.pdb")
    working_pose.dump_pdb(output_pdb)
    print(f"Saved docked pose to: {output_pdb}")
    total = dockscore(working_pose)
    score_data = {
        "run_id": N,
        "output_name": os.path.basename(output_pdb),
        "total_score": total,
        "score_per_res": total / working_pose.size()
    }
    emap = working_pose.energies().total_energies()
    for k in emap:
        score_data[k.name()] = emap[k]

    all_scores.append(score_data)
    print(f" Run {N} finished in {round(time.time() - t0, 2)} seconds.")

df = pd.DataFrame(all_scores)
df.to_csv(scorefile, index=False)
print(f" All scores + energy breakdown saved to: {scorefile}")
'''
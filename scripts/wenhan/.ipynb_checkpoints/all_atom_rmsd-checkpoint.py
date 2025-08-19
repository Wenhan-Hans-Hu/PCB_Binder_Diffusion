#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 11:15:49 2025

@author: Wenhan Hu
"""
'''
Example Use:
    python ligand_rmsd.py \
    --pdb 1.pdb \
    --ref_pdb 2.pdb \
    --params LIG.params \
    --ligand LIG
'''
import os
import argparse
import csv
import pyrosetta
from pyrosetta import pose_from_file, Pose
from pyrosetta.rosetta.core.scoring import all_atom_rmsd
from pyrosetta.rosetta.protocols.toolbox import CA_superimpose
from pyrosetta.rosetta.std import list_unsigned_long_t
from pyrosetta.rosetta.core.scoring import CA_rmsd

# Path to CSV
csv_path = os.environ.get("RMSD_CSV_PATH", "all_rmsd_results.csv")

def main():
    parser = argparse.ArgumentParser(description="Accurate ligand RMSD after protein alignment.")
    parser.add_argument("--target_pdb", type=str, required=True)
    parser.add_argument("--ref_pdb", type=str, required=True)
    parser.add_argument("--params", nargs="+", required=True)
    args = parser.parse_args()

    # Init PyRosetta
    extra_res = ' '.join(f"-extra_res_fa {p}" for p in args.params)
    pyrosetta.init(f"{extra_res} -run:preserve_header")
    
    # Load poses reference, but i won't call them directly
    pose_target = pose_from_file(args.target_pdb)
    pose_ref = pose_from_file(args.ref_pdb)

    # ###Aligned by ligand part (1,2)
    # pose_1 = pose_target.clone()
    # lig_1 = pose_1.size()
    # assert pose_1.residue(lig_1).is_ligand()   
    # ligand_pose_1 = Pose()
    # ligand_pose_1.append_residue_by_bond(pose_1.residue(lig_1))
    
    # pose_2 = pose_ref.clone()
    # lig_2 = pose_2.size()
    # assert pose_2.residue(lig_2).is_ligand() 
    # ligand_pose_2 = Pose()
    # ligand_pose_2.append_residue_by_bond(pose_2.residue(lig_2))
    
    # # Align ligands ( This line might be wrong )
    # before = CA_rmsd(ligand_pose_1, ligand_pose_2)
    # print(f"Ligand Cα RMSD before alignment: {before:.4f} Å")
    
    # CA_superimpose(ligand_pose_1, ligand_pose_2)
    # ligand_only_rmsd = CA_rmsd(ligand_pose_1, ligand_pose_2)
    # print(f"Ligand Cα RMSD after alignment: {ligand_only_rmsd:.4f} Å")

    # res_list = list_unsigned_long_t()
    # res_list.append(lig1)
    # ligand_only_rmsd = all_atom_rmsd(pose2, pose1, res_list)

    ### Aligned by protein part (a,b)
    pose_a = pose_target.clone()
    lig_a = pose_a.size()
    assert pose_a.residue(lig_a).is_ligand()
    ligand_pose_a = Pose()
    ligand_pose_a.append_residue_by_bond(pose_a.residue(lig_a))
    
    pose_b = pose_ref.clone()
    lig_b = pose_b.size()
    assert pose_b.residue(lig_b).is_ligand()
    ligand_pose_b = Pose()
    ligand_pose_b.append_residue_by_bond(pose_b.residue(lig_b))

    Before = CA_rmsd(ligand_pose_a, ligand_pose_b)
    print(f"Protein Cα RMSD before alignment: {Before:.4f} Å")    
    CA_superimpose(pose_a, pose_b)
    After = CA_rmsd(ligand_pose_a, ligand_pose_b)
    print(f"Protein Cα RMSD after alignment: {After:.4f} Å")
    
    # RMSD between aligned ligand poses
    rmsd = all_atom_rmsd(ligand_pose_a, ligand_pose_b)
    # print(f"Ligand-only RMSD (before alignment): {ligand_only_rmsd:.4f} Å")
    print(f"Ligand all-atom RMSD (aligned proteins): {rmsd:.4f} Å")
    

    # Output CSV
    # write_header = not os.path.exists(csv_path)
    write_header = not os.path.isfile(csv_path) or os.stat(csv_path).st_size == 0

    with open(csv_path, "a", newline="") as csvfile:
        writer = csv.writer(csvfile)
        if write_header:
            writer.writerow(["Target", "Reference", "Ligand_ResID", "RMSD_Protein_Aligned"])
        writer.writerow([
            os.path.basename(args.target_pdb),
            os.path.basename(args.ref_pdb),
            lig_a,
            f"{rmsd:.4f}",
        ])


if __name__ == "__main__":
    main()

'''
def main():
    parser = argparse.ArgumentParser(description="Accurate ligand RMSD after protein alignment.")
    parser.add_argument("--target_pdb", type=str, required=True)
    parser.add_argument("--ref_pdb", type=str, required=True)
    parser.add_argument("--params", nargs="+", required=True)
    args = parser.parse_args()

    # Init PyRosetta
    extra_res = ' '.join(f"-extra_res_fa {p}" for p in args.params)
    pyrosetta.init(f"{extra_res} -run:preserve_header")
    
    ### Un aligned by protein part
    target = pose_from_file(args.target_pdb)
    ref = pose_from_file(args.ref_pdb)
    lig1 = target.size()
    lig2 = ref.size()
    assert target.residue(lig1).is_ligand()
    assert ref.residue(lig2).is_ligand()    
    res_list = list_unsigned_long_t()
    res_list.append(lig1)
    ligand_only_rmsd = all_atom_rmsd(target, ref, res_list)

    ### Aligned by protein part
    # Load full poses
    pose_target = pose_from_file(args.target_pdb)
    pose_ref = pose_from_file(args.ref_pdb)

    # Ligand assumed to be last residue
    lig_resnum_target = pose_target.size()
    lig_resnum_ref = pose_ref.size()

    assert pose_target.residue(lig_resnum_target).is_ligand()
    assert pose_ref.residue(lig_resnum_ref).is_ligand()

    # Align proteins
    CA_superimpose(pose_target, pose_ref)

    # Extract ligand-only poses (1 residue each)
    ligand_pose_target = Pose()
    ligand_pose_target.append_residue_by_bond(pose_target.residue(lig_resnum_target))

    ligand_pose_ref = Pose()
    ligand_pose_ref.append_residue_by_bond(pose_ref.residue(lig_resnum_ref))

    # RMSD between aligned ligand poses
    rmsd = all_atom_rmsd(ligand_pose_target, ligand_pose_ref)
    print(f"Ligand-only RMSD (before alignment): {ligand_only_rmsd:.4f} Å")
    print(f"Ligand all-atom RMSD (aligned proteins): {rmsd:.4f} Å")
    

    # Output CSV
    write_header = not os.path.exists(csv_path)
    with open(csv_path, "a", newline="") as csvfile:
        writer = csv.writer(csvfile)
        if write_header:
            writer.writerow(["Target", "Reference", "Ligand_ResID", "RMSD_Protein_Aligned", "RMSD_Ligand_Only"])
        writer.writerow([
            os.path.basename(args.target_pdb),
            os.path.basename(args.ref_pdb),
            lig_resnum_target,
            f"{rmsd:.4f}",
            f"{ligand_only_rmsd:.4f}"
        ])
'''
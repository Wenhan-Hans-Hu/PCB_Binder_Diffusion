#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 18:01:02 2025

@author: wh621
"""

import argparse
import os
import csv
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS

# Parse input arguments
parser = argparse.ArgumentParser(description="RDKit RMSD calculator")
parser.add_argument("--ref_pdb", required=True, help="Reference PDB file")
parser.add_argument("--target_pdb", required=True, help="Target PDB file")
parser.add_argument("--params", nargs="+", type=str, help=".params file(s) for ligand/cofactors")
args = parser.parse_args()

ref_pdb = args.ref_pdb
target_pdb = args.target_pdb
csv_path = os.environ.get("RMSD_CSV_PATH")

'''
mol_tar = Chem.MolFromPDBFile(tar.pdb, removeHs=True)
mol_ref = Chem.MolFromPDBFile("ref.pdb", removeHs=False)
if mol_tar is None or mol_ref is None:
    raise ValueError("Failed to load one or both PDB files.")
rmsd = AllChem.GetBestRMS(mol_tar, mol_ref)
print(f"RMSD between tar.pdb and ref.pdb: {rmsd:.3f} Å")
'''

# Load molecules
mol_tar = Chem.MolFromPDBFile(target_pdb, removeHs=True)
mol_ref = Chem.MolFromPDBFile(ref_pdb, removeHs=True)

if mol_ref is None or mol_tar is None:
    raise RuntimeError(f"Failed to load one or both PDB files: {ref_pdb}, {target_pdb}")

# # Try to map atoms using Maximum Common Substructure
# mcs_result = rdFMCS.FindMCS([mol_ref, mol_tar])
# mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
# match_ref = mol_ref.GetSubstructMatch(mcs_mol)
# match_tar = mol_tar.GetSubstructMatch(mcs_mol)

# if not match_ref or not match_tar:
#     raise RuntimeError(f"Failed to find MCS between {ref_pdb} and {target_pdb}")

# atom_map = list(zip(match_tar, match_ref))
rmsd = AllChem.GetBestRMS(mol_tar, mol_ref)

# Output
ref_name = os.path.basename(ref_pdb)
tar_name = os.path.basename(target_pdb)

print(f"RMSD between {ref_name} and {tar_name}: {rmsd:.3f} Å")

# Append to CSV
header = ["ref_pdb", "target_pdb", "rmsd"]
row = [ref_name, tar_name, f"{rmsd:.3f}"]

file_exists = os.path.isfile(csv_path)
with open(csv_path, "a", newline="") as f:
    writer = csv.writer(f)
    if not file_exists:
        writer.writerow(header)
    writer.writerow(row)
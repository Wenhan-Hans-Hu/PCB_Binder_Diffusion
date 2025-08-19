#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFreeSASA
from rdkit.Chem import AllChem

if len(sys.argv) < 2:
    print("Usage: python sasa_batch_single.py <pdb_file>")
    sys.exit(1)

input_pdb = sys.argv[1]
input_pdb_name = os.path.splitext(os.path.basename(input_pdb))[0]
ligand_name = "CYC"  # You can change this as needed
output_file = "ligand_sasa_output.xlsx"

def splitProteinLigandFromPDB_All(pdb_path, ligand_resname="res"):
    mol = Chem.rdmolfiles.MolFromPDBFile(pdb_path, removeHs=False, sanitize=False)
    if mol is None:
        raise ValueError("Failed to read PDB file.")

    prot_atoms = []
    lig_atom_groups = {}

    for atom in mol.GetAtoms():
        resinfo = atom.GetPDBResidueInfo()
        if resinfo is None:
            continue
        resname = resinfo.GetResidueName().strip()
        chain_id = resinfo.GetChainId().strip()
        res_id = resinfo.GetResidueNumber()
        key = (chain_id, res_id)

        if resname == ligand_resname:
            if key not in lig_atom_groups:
                lig_atom_groups[key] = []
            lig_atom_groups[key].append(atom.GetIdx())
        else:
            prot_atoms.append(atom.GetIdx())

    if len(prot_atoms) == 0:
        raise ValueError("No protein atoms found.")
    if len(lig_atom_groups) == 0:
        raise ValueError(f"No ligand atoms found for residue '{ligand_resname}'.")

    prot = Chem.PathToSubmol(mol, prot_atoms)
    lig_list = [Chem.PathToSubmol(mol, atom_idxs) for atom_idxs in lig_atom_groups.values()]
    return prot, lig_list, list(lig_atom_groups.keys())

def calcSASA(mol):
    ptable = Chem.GetPeriodicTable()
    radii = [ptable.GetRvdw(atom.GetAtomicNum()) for atom in mol.GetAtoms()]
    rdFreeSASA.CalcSASA(mol, radii)
    return mol

def calcEachLigandSASA(complex_pdbfile, ligand_resname="res"):
    prot, ligands, ligand_keys = splitProteinLigandFromPDB_All(complex_pdbfile, ligand_resname)

    results = []
    for i, (lig, key) in enumerate(zip(ligands, ligand_keys)):
        lig = calcSASA(lig)
        lig_sasa = sum(float(a.GetProp("SASA")) for a in lig.GetAtoms())

        comp = Chem.CombineMols(prot, lig)
        comp = Chem.AddHs(comp, addCoords=True)
        comp = calcSASA(comp)

        comp_lig = Chem.GetMolFrags(comp, asMols=True, sanitizeFrags=False)[-1]
        bound_sasa = sum(float(a.GetProp("SASA")) for a in comp_lig.GetAtoms())

        results.append({
            "file_name": os.path.basename(complex_pdbfile),
            "ligand_no": i + 1,
            "chain": key[0],
            "Free ligand SASA": lig_sasa,
            "Bound ligand SASA": bound_sasa,
            "Ratio": bound_sasa / lig_sasa if lig_sasa != 0 else 0.0
        })

    return results

# Load or create the output DataFrame
if os.path.exists(output_file):
    sasa_df = pd.read_excel(output_file)
else:
    sasa_df = pd.DataFrame(columns=["file_name", "ligand_no", "chain", "Free ligand SASA", "Bound ligand SASA", "Ratio"])

# Calculate and append new results
try:
    results = calcEachLigandSASA(input_pdb, ligand_name)
    sasa_df = pd.concat([sasa_df, pd.DataFrame(results)], ignore_index=True)
    sasa_df.to_excel(output_file, index=False)
    print(f"Results saved to {output_file}")
except Exception as e:
    print(f"Error: {e}")


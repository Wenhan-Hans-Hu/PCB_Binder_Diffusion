#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 19:23:16 2023

@author: ikalvet
"""
import pyrosetta as pyr
import pyrosetta.rosetta
import os, sys
import pandas as pd

SCRIPT_PATH = os.path.dirname(__file__)
sys.path.append(f"{SCRIPT_PATH}/../../utils")
import no_ligand_repack
import scoring_utils
import design_utils


comparisons = {'<=': '__le__',
               '<': '__lt__',
               '>': '__gt__',
               '>=': '__ge__',
               '=': '__eq__'}


def fix_scorefxn(sfxn, allow_double_bb=False):
    opts = sfxn.energy_method_options()
    opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    opts.hbond_options().bb_donor_acceptor_check(not allow_double_bb)
    sfxn.set_energy_method_options(opts)


def score_design(pose, sfx, catres):
    df_scores = pd.DataFrame()

    # Adding Rosetta scores to df
    sfx(pose)
    for k in pose.scores:
        df_scores.at[0, k] = pose.scores[k]

    # Calculating CST score
    if pose.constraint_set().has_constraints():
        sfx.set_weight(pyrosetta.rosetta.core.scoring.score_type_from_name("atom_pair_constraint"), 1.5)
        sfx.set_weight(pyrosetta.rosetta.core.scoring.score_type_from_name("angle_constraint"), 1.0)
        sfx.set_weight(pyrosetta.rosetta.core.scoring.score_type_from_name("dihedral_constraint"), 1.0)
        sfx(pose)
        df_scores.at[0, "all_cst"] = sum([pose.scores[s] for s in pose.scores if "constraint" in s])

    # Calculating ddg
    # LUKAS: Absolutely not! This measures change in ligand binding upon
    # mutation of the catalytic residue, which is unnecessary when not doing
    # heme. This whole section has thus been commented into non-existence.
    
#    from scoring_utils import calculate_ddg

    fix_scorefxn(sfx)

#    df_scores.at[0, 'corrected_ddg'] = calculate_ddg(pose, sfx)

    # Calculating relative ligand SASA
    # First figuring out what is the path to the ligand PDB file
    ligand_seqpos = pose.size()
    assert pose.residue(ligand_seqpos).is_ligand()

    ligand_pose = pyrosetta.rosetta.core.pose.Pose()
    pyrosetta.rosetta.core.pose.append_subpose_to_pose(ligand_pose, pose, pose.size(), pose.size(), 1)

    free_ligand_sasa = scoring_utils.getSASA(ligand_pose, resno=1)
    ligand_sasa = scoring_utils.getSASA(pose, resno=ligand_seqpos)
    df_scores.at[0, 'L_SASA'] = ligand_sasa / free_ligand_sasa


    # Using a custom function to find HBond partners of the COO groups. LUKAS: Heme only!
#    for n in range(1, 5):
#        df_scores.at[0, f"O{n}_hbond"] = scoring_utils.find_hbonds_to_residue_atom(pose, ligand_seqpos, f"O{n}")


    # Checking if both COO groups on heme are hbonded. LUKAS: Heme only!!!
#    if any([df_scores.at[0, x] > 0.0 for x in ['O1_hbond', 'O3_hbond']]) and any([df_scores.at[0, x] > 0.0 for x in ['O2_hbond', 'O4_hbond']]):
#        df_scores.at[0, 'COO_hbond'] = 1.0
#    else:
#        df_scores.at[0, 'COO_hbond'] = 0.0

    # Calculating ContactMolecularSurface
    cms = pyrosetta.rosetta.protocols.simple_filters.ContactMolecularSurfaceFilter()
    lig_sel = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(ligand_seqpos)
    protein_sel = pyrosetta.rosetta.core.select.residue_selector.ChainSelector("A")
    cms.use_rosetta_radii(True)
    cms.distance_weight(0.5)
    cms.selector1(protein_sel)
    cms.selector2(lig_sel)
    df_scores.at[0, "cms"] =  cms.compute(pose)  # ContactMolecularSurface
    df_scores.at[0, "cms_per_atom"] =  df_scores.at[0, "cms"] / pose.residue(ligand_seqpos).natoms()  # ContactMolecularSurface per atom

    ## Calculating shape complementarity
    sc = pyrosetta.rosetta.protocols.simple_filters.ShapeComplementarityFilter()
    sc.use_rosetta_radii(True)
    sc.selector1(protein_sel)
    sc.selector2(lig_sel)
    df_scores.at[0, "sc"] = sc.score(pose)

    # Running no-ligand-repack. LUKAS: Knocked out this section due to bugs being thrown during Fe binder design. Possibly needs uncommenting later. Fe cannot renerate catalytic resnos so this will always fail. LUKAS: I think I'm wrong. Ligand_resno comes from ligand_seqpos comes from pose.size().
    nlr_scores = no_ligand_repack.no_ligand_repack(pose, pyr.get_fa_scorefxn(), ligand_resno=ligand_seqpos)
    for k in nlr_scores.keys():
        df_scores.at[0, k] = nlr_scores.iloc[0][k]

    df_scores.at[0, "score_per_res"] = df_scores.at[0, "total_score"]/pose.size()

    return df_scores


def filter_scores(scores):
    """
    Filters are defined in this importable module
    """
    filtered_scores = scores.copy()

    for s in filters.keys():
        if filters[s] is not None and s in scores.keys():
            val = filters[s][0]
            sign = comparisons[filters[s][1]]
            filtered_scores =\
              filtered_scores.loc[(filtered_scores[s].__getattribute__(sign)(val))]
            n_passed = len(scores.loc[(scores[s].__getattribute__(sign)(val))])
            print(f"{s:<24} {filters[s][1]:<2} {val:>7.3f}: {len(filtered_scores)} "
                  f"designs left. {n_passed} pass ({(n_passed/len(scores))*100:.0f}%).")
    return filtered_scores


filters = {#"all_cst": [1.5, "<="],
           "L_SASA": [0.2, "<="],
           #"COO_hbond": [1.0, "<="],
           "cms_per_atom": [5.0, ">="],
           #"corrected_ddg": [50.0, "<="],
           "nlr_totrms": [0.8, "<="]}
           #"nlr_SR1_rms": [0.6, "<="]}

# LUKAS: I certainly hope this variable isn't being called, but in case it is, I've dropped a reference to "N4" because TNT does not have 4 atoms. 
# UPDATE: This variable is being called. It throws an error if you remove it because it needs it, seemingly, to determine what residues are in the ligand pocket and which are too far from the ligand to be worth redesigning. If you delete the variable the script will throw an AssertionError, if you have no atoms in it it wil continuously say "NAN occurred in hbonding calculations." This means that you need to specify a set of atoms in your ligand with which to do the calculations with. For TNT, I use Nitrogen atoms 1-3, but you need to choose a set actually exists in your ligand, or it will throw "atom x does not exist in ligand!"
# Yes, this is very stupid and bad programming. Life goes on.
align_atoms = ["N1", "N2", "N3", "N4"]

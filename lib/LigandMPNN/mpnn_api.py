import argparse
import sys
import copy
import io
import torch
import random
import json
import numpy as np
import pandas as pd
import os.path
import data_utils
from data_utils import get_seq_rec, get_score, parse_PDB, featurize, write_full_PDB
from model_utils import ProteinMPNN
import time
SCRIPT_PATH = os.path.dirname(__file__)


restype_1to3 = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL', 'X': 'UNK'}
restype_STRtoINT = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19, 'X': 20}
restype_INTtoSTR = {0: 'A', 1: 'C', 2: 'D', 3: 'E', 4: 'F', 5: 'G', 6: 'H', 7: 'I', 8: 'K', 9: 'L', 10: 'M', 11: 'N', 12: 'P', 13: 'Q', 14: 'R', 15: 'S', 16: 'T', 17: 'V', 18: 'W', 19: 'Y', 20: 'X'}
alphabet = list(restype_STRtoINT)

class MPNNRunner(object):
    def __init__(self, model_type, checkpoint_path=None, ligand_mpnn_use_side_chain_context=False, seed=None, verbose=False):
        """
            model_type (str) :: must be one of:
                ['protein_mpnn', 'ligand_mpnn', 'per_residue_label_membrane_mpnn',
                 'global_label_membrane_mpnn', 'soluble_mpnn']
            checkpoint_path (str)
            ligand_mpnn_use_side_chain_context (bool)
            seed (int)
            verbose (bool)
        """

        #fix seeds
        if seed is not None:
            self.seed=seed
        else:
            self.seed=int(np.random.randint(0, high=99999, size=1, dtype=int)[0])
        
        torch.manual_seed(self.seed)
        random.seed(self.seed)
        np.random.seed(self.seed)

        __checkpoints = {"protein_mpnn": f"{SCRIPT_PATH}/model_params/proteinmpnn_v_48_020.pt",
                         "ligand_mpnn": f"{SCRIPT_PATH}/model_params/ligandmpnn_v_32_010_25.pt",
                         "per_residue_label_membrane_mpnn": f"{SCRIPT_PATH}/model_params/per_residue_label_membrane_mpnn_v_48_020.pt",
                         "global_label_membrane_mpnn": f"{SCRIPT_PATH}/model_params/global_label_membrane_mpnn_v_48_020.pt",
                         "soluble_mpnn": f"{SCRIPT_PATH}/model_params/solublempnn_v_48_020.pt"}

        if checkpoint_path is None:
            assert model_type in __checkpoints.keys(), "invalid model_type input"
            self.__checkpoint_path = __checkpoints[model_type]
        else:
            print(f"User provided checkpoint: {checkpoint_path}.")
            print("Warning! No checks are done to make sure that this is the right checkpoint for your task.")
            self.__checkpoint_path = checkpoint_path

        self.__model_type = model_type
        self.ligand_mpnn_use_side_chain_context = ligand_mpnn_use_side_chain_context

        print(f"Using {model_type} model from checkpoint: {self.__checkpoint_path}")

        self.device = torch.device("cuda:0" if (torch.cuda.is_available()) else "cpu")

        checkpoint = torch.load(self.__checkpoint_path, map_location=self.device)

        if model_type == "ligand_mpnn":
            self.atom_context_num = checkpoint["atom_context_num"]
            k_neighbors=32
        else:
            self.atom_context_num = 1
            self.ligand_mpnn_use_side_chain_context = False
            k_neighbors=checkpoint["num_edges"]

        self.model = ProteinMPNN(node_features=128,
                        edge_features=128,
                        hidden_dim=128,
                        num_encoder_layers=3,
                        num_decoder_layers=3,
                        k_neighbors=k_neighbors,
                        device=self.device,
                        atom_context_num=self.atom_context_num,
                        model_type=model_type,
                        ligand_mpnn_use_side_chain_context=self.ligand_mpnn_use_side_chain_context)

        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.model.to(self.device)
        self.model.eval()

        self.verbose = verbose
        pass

    def checkpoint(self):
        return self.__checkpoint

    @property
    def model_type(self):
        return self.__model_type

    class MPNN_Input(object):
        def __init__(self, obj=None):
            ## This just crudely replicates the args object with its attribute namespace
            ## TODO: implement taking another object of this class as input, and copy the attribute values
            self.structure = None
            self.omit_AA = []
            self.tied_positions = None
            self.bias_by_res = None
            self.bias_AA = None
            self.temperature = None
            self.name = "pose_0000"
            self.max_length = 20000
            self.pdb = None  # string representation of a PDB
            self.fixed_residues = []
            self.redesigned_residues = []
            self.parse_these_chains_only = None
            self.bias_AA_per_residue = None
            self.omit_AA_per_residue = None
            self.transmembrane_buried = None
            self.transmembrane_interface = None
            self.global_transmembrane_label = None
            self.chains_to_design = None
            self.symmetry_residues = None
            self.symmetry_weights = None
            self.homo_oligomer = None
            self.verbose = None
            self.ligand_mpnn_cutoff_for_score = 8.0
            self.ligand_mpnn_use_atom_context = True
            self.parse_atoms_with_zero_occupancy = False
            self.batch_size = None
            self.number_of_batches = None

            self.zero_indexed = None
            self.force_hetatm = None

            # Cloning the values of input object, if given
            if obj is not None:
                assert isinstance(obj, MPNNRunner.MPNN_Input)
                for attr in obj.__dir__():
                    if attr[:2] == "__":
                        continue
                    self.__setattr__(attr, copy.deepcopy(obj.__getattribute__(attr)))
            pass

        def copy(self):
            return copy.deepcopy(self)


    def run(self, input_obj, **kwargs):

        if self.verbose:
            print("Designing this PDB:", input_obj.name)
        fixed_residues = input_obj.fixed_residues
        redesigned_residues = input_obj.redesigned_residues

        bias_AA = torch.zeros([21], device=self.device, dtype=torch.float32)
        if input_obj.bias_AA:
            assert isinstance(input_obj.bias_AA, dict)
            for AA, bias in input_obj.bias_AA.items():
                bias_AA[restype_STRtoINT[AA]] = bias

        #make array to indicate which amino acids need to be omitted [21]
        omit_AA = torch.tensor(np.array([AA in input_obj.omit_AA for AA in alphabet]).astype(np.float32), device=self.device)

        #parse PDB file
        protein_dict, backbone, other_atoms, icodes, water_atoms = parse_PDB(io.StringIO(input_obj.pdb),
                                                                            device=self.device, 
                                                                            chains=input_obj.parse_these_chains_only,
                                                                            parse_all_atoms=self.ligand_mpnn_use_side_chain_context,
                                                                            parse_atoms_with_zero_occupancy=input_obj.parse_atoms_with_zero_occupancy)

        #----
        #make chain_letter + residue_idx + insertion_code mapping to integers
        R_idx_list = list(protein_dict["R_idx"].cpu().numpy())
        chain_letters_list = list(protein_dict["chain_letters"])
        encoded_residues = []
        for i in range(len(R_idx_list)):
            tmp = str(chain_letters_list[i]) + str(R_idx_list[i]) + icodes[i]
            encoded_residues.append(tmp)
        encoded_residue_dict = dict(zip(encoded_residues, range(len(encoded_residues))))
        encoded_residue_dict_rev = dict(zip(list(range(len(encoded_residues))), encoded_residues))

        bias_AA_per_residue = torch.zeros([len(encoded_residues),21], device=self.device, dtype=torch.float32)
        if input_obj.bias_AA_per_residue:    
            bias_dict = input_obj.bias_AA_per_residue
            for residue_name, v1 in bias_dict.items():
                if residue_name in encoded_residues:
                    i1 = encoded_residue_dict[residue_name]
                    for amino_acid, v2 in v1.items():
                        if amino_acid in alphabet:
                            j1 = restype_STRtoINT[amino_acid]
                            bias_AA_per_residue[i1,j1] = v2
        #----

        omit_AA_per_residue = torch.zeros([len(encoded_residues),21], device=self.device, dtype=torch.float32)
        if input_obj.omit_AA_per_residue:    
            omit_dict = input_obj.omit_AA_per_residue
            for residue_name, v1 in omit_dict.items():
                if residue_name in encoded_residues:
                    i1 = encoded_residue_dict[residue_name]
                    for amino_acid in v1:
                        if amino_acid in alphabet:
                            j1 = restype_STRtoINT[amino_acid]
                            omit_AA_per_residue[i1,j1] = 1.0
        #----


        fixed_positions = torch.tensor([int(item not in fixed_residues) for item in encoded_residues], device=self.device)
        redesigned_positions = torch.tensor([int(item not in redesigned_residues) for item in encoded_residues], device=self.device)
        #----

        #specify which residues are buried for checkpoint_per_residue_label_membrane_mpnn model
        if input_obj.transmembrane_buried:
            buried_residues = input_obj.transmembrane_buried
            buried_positions = torch.tensor([int(item in buried_residues) for item in encoded_residues], device=self.device)
        else:
            buried_positions = torch.zeros_like(fixed_positions)
        #----

        if input_obj.transmembrane_interface:
            interface_residues = input_obj.transmembrane_interface
            interface_positions = torch.tensor([int(item in interface_residues) for item in encoded_residues], device=self.device)
        else:
            interface_positions = torch.zeros_like(fixed_positions)
        #----
        protein_dict["membrane_per_residue_labels"] = 2*buried_positions*(1-interface_positions) + 1*interface_positions*(1-buried_positions)

        if self.model_type == "global_label_membrane_mpnn":
            protein_dict["membrane_per_residue_labels"] = input_obj.global_transmembrane_label + 0*fixed_positions

        #specify which chains need to be redesigned
        if isinstance(input_obj.chains_to_design, str):
            chains_to_design_list = input_obj.chains_to_design.split(",")
        elif isinstance(input_obj.chains_to_design, list):
            chains_to_design_list = input_obj.chains_to_design
        else:
            chains_to_design_list = protein_dict["chain_letters"]
        chain_mask = torch.tensor(np.array([item in chains_to_design_list for item in protein_dict["chain_letters"]],dtype=np.int32), device=self.device)
        #----

        #create chain_mask to notify which residues are fixed (0) and which need to be designed (1)
        if redesigned_residues:
            protein_dict["chain_mask"] = chain_mask*(1-redesigned_positions)
        elif fixed_residues:
            protein_dict["chain_mask"] = chain_mask*fixed_positions
        else:
            protein_dict["chain_mask"] = chain_mask

        if self.verbose:
            PDB_residues_to_be_redesigned = [encoded_residue_dict_rev[item] for item in range(protein_dict["chain_mask"].shape[0]) if protein_dict["chain_mask"][item]==1]
            PDB_residues_to_be_fixed = [encoded_residue_dict_rev[item] for item in range(protein_dict["chain_mask"].shape[0]) if protein_dict["chain_mask"][item]==0]
            print("These residues will be redesigned: ", PDB_residues_to_be_redesigned)
            print("These residues will be fixed: ", PDB_residues_to_be_fixed)
        #----


        #----
        #specify which residues are linked
        if input_obj.symmetry_residues:
            symmetry_residues_list_of_lists = [x.split(',') for x in input_obj.symmetry_residues.split('|')]
            remapped_symmetry_residues=[]
            for t_list in symmetry_residues_list_of_lists:
                tmp_list=[]
                for t in t_list:
                    tmp_list.append(encoded_residue_dict[t])
                remapped_symmetry_residues.append(tmp_list) 
        else:
            remapped_symmetry_residues=[[]]
        #----

        #specify linking weights
        if input_obj.symmetry_weights:
            symmetry_weights = [[float(item) for item in x.split(',')] for x in input_obj.symmetry_weights.split('|')]
        else:
            symmetry_weights = [[]]
        #----

        if input_obj.homo_oligomer:
            if input_obj.verbose:
                print("Designing HOMO-OLIGOMER")
            chain_letters_set = list(set(chain_letters_list))
            reference_chain = chain_letters_set[0]
            lc = len(reference_chain)
            residue_indices = [item[lc:] for item in encoded_residues if item[:lc]==reference_chain]
            remapped_symmetry_residues=[]
            symmetry_weights = []
            for res in residue_indices:
                tmp_list=[]
                tmp_w_list=[]
                for chain in chain_letters_set:
                    name = chain+res
                    tmp_list.append(encoded_residue_dict[name])
                    tmp_w_list.append(1.0)
                remapped_symmetry_residues.append(tmp_list)
                symmetry_weights.append(tmp_w_list)

        out_dict = {}

        with torch.no_grad():
            #run featurize to remap R_idx and add batch dimension
            feature_dict = featurize(protein_dict,
                                    cutoff_for_score=input_obj.ligand_mpnn_cutoff_for_score, 
                                    use_atom_context=input_obj.ligand_mpnn_use_atom_context,
                                    number_of_ligand_atoms=self.atom_context_num,
                                    model_type=self.model_type)
            feature_dict["batch_size"] = input_obj.batch_size
            B, L, _, _ = feature_dict["X"].shape #batch size should be 1 for now.
            #----

            #add additional keys to the feature dictionary
            feature_dict["temperature"] = input_obj.temperature
            feature_dict["bias"] = (-1e8*omit_AA[None,None,:]+bias_AA).repeat([1,L,1])+bias_AA_per_residue[None]-1e8*omit_AA_per_residue[None]
            feature_dict["symmetry_residues"] = remapped_symmetry_residues
            feature_dict["symmetry_weights"] = symmetry_weights
            #----

            sampling_probs_list = []
            log_probs_list = []
            decoding_order_list = []
            S_list = []
            loss_list = []
            loss_per_residue_list = []
            loss_XY_list = []
            for _ in range(input_obj.number_of_batches):
                feature_dict["randn"] = torch.randn([feature_dict["batch_size"], feature_dict["mask"].shape[1]], device=self.device)
                #main step-----
                output_dict = self.model.sample(feature_dict)

                #compute confidence scores
                loss, loss_per_residue = get_score(output_dict["S"], output_dict["log_probs"], feature_dict["mask"]*feature_dict["chain_mask"])
                if self.model_type == "ligand_mpnn":
                    combined_mask = feature_dict["mask"]*feature_dict["mask_XY"]*feature_dict["chain_mask"]
                else:
                    combined_mask = feature_dict["mask"]*feature_dict["chain_mask"]
                loss_XY, _ = get_score(output_dict["S"], output_dict["log_probs"], combined_mask)
                #-----
                S_list.append(output_dict["S"])
                log_probs_list.append(output_dict["log_probs"])
                sampling_probs_list.append(output_dict["sampling_probs"])
                decoding_order_list.append(output_dict["decoding_order"])
                loss_list.append(loss)
                loss_per_residue_list.append(loss_per_residue)
                loss_XY_list.append(loss_XY)
            S_stack = torch.cat(S_list, 0)
            log_probs_stack = torch.cat(log_probs_list, 0)
            sampling_probs_stack = torch.cat(sampling_probs_list, 0)
            decoding_order_stack = torch.cat(decoding_order_list, 0)
            loss_stack = torch.cat(loss_list, 0)
            loss_per_residue_stack = torch.cat(loss_per_residue_list, 0)
            loss_XY_stack = torch.cat(loss_XY_list, 0)
            rec_mask = feature_dict["mask"][:1]*feature_dict["chain_mask"][:1]
            rec_stack = get_seq_rec(feature_dict["S"][:1], S_stack, rec_mask)

            
            #make input sequence string separated by / between different chains
            native_seq = "".join([restype_INTtoSTR[AA] for AA in feature_dict["S"][0].cpu().numpy()])
            seq_np = np.array(list(native_seq))
            seq_out_str = []
            for mask in protein_dict['mask_c']:
                seq_out_str += list(seq_np[mask.cpu().numpy()])
                seq_out_str += ['/']
            seq_out_str = "".join(seq_out_str)[:-1]
            native_seq_str = seq_out_str
            #------

            out_dict["generated_sequences_int"] = S_stack.detach().cpu().numpy()
            out_dict["generated_sequences"] = ["".join([restype_INTtoSTR[AA] for AA in s_ix]) for s_ix in S_stack.detach().cpu().numpy()]
            out_dict["scores"] = [np.exp(-L_ix) for L_ix in loss_stack.detach().cpu().numpy()]
            out_dict["sampling_probs"] = sampling_probs_stack.detach().cpu().numpy()
            out_dict["sampling_probs_dict"] = [[{restype_INTtoSTR[i]: v for i,v in enumerate(res_probs)} for res_probs in seqprobs] for seqprobs in out_dict["sampling_probs"]]
            out_dict["log_probs"] = log_probs_stack.detach().cpu().numpy()
            out_dict["decoding_order"] = decoding_order_stack.detach().cpu().numpy()
            out_dict["native_sequence"] = feature_dict["S"][0].detach().cpu().numpy()
            out_dict["mask"] = feature_dict["mask"][0].detach().cpu().numpy()
            out_dict["chain_mask"] = feature_dict["chain_mask"][0].detach().cpu().numpy()
            out_dict["seed"] = self.seed
            out_dict["temperature"] = input_obj.temperature

            return out_dict

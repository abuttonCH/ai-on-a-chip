import numpy as np
from rdkit.Chem import MolToSmiles as MS
from rdkit.Chem import MolFromSmiles as SM
import os, sys
import multiprocessing as mp
import time
import argparse

class retrieve_building_blocks:
    def __init__(self,building_block_set,output_path,limit):
        self.building_block_set = np.asarray([b.decode("UTF-8") for b in building_block_set])
        self.output_path = output_path
        #limit the number of molecules that are read in
        self.limit = limit

    def search_for_bb_string_compare(self,reactant_mol):
        #find id of matching smiles
        selected_indices = np.where(self.building_block_set == reactant_mol)[0]
        if len(selected_indices) != 0:
            selected_bb = self.building_block_set[selected_indices[0]]
        else:
            selected_bb = None
        return selected_bb

    def retrieve_bb_pairs(self,decomp_data):
        output_file = open(self.output_path,"w+")
        output_file.write("ID|product molecule|decomposition reaction|retrieved reactants\n")

        for index,decomp_entry in enumerate(decomp_data):
            #limit the number of molecules that are read in
            if self.limit != None:
                if index > int(self.limit):
                    break
            product_id = decomp_entry[0]
            LSTM_product = decomp_entry[1]
            reaction_name = decomp_entry[2]
            reactant_pair = decomp_entry[-1].split(".")

            retrieved_building_blocks = []
            for mol in reactant_pair:
                #perform retrieval
                select_bb = self.search_for_bb_string_compare(mol)
                if select_bb != None:
                    retrieved_building_blocks.append(select_bb)
            #only logs results if all reactants were retrieved
            if len(retrieved_building_blocks) == len(reactant_pair):
                output_file.write("|".join(decomp_entry)+"\n")
            else:
                pass
        output_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--decomp')
    parser.add_argument('--mol_db')
    parser.add_argument('--out')
    parser.add_argument('--limit')
    args = parser.parse_args()

    #inputs
    decompo_set_file = open(args.decomp,"r") #decomposed molecules
    mol_set = np.load(args.mol_db)[:,1] #molecular database
    mol_limit = args.limit #limit on the number of molecules
    output_path = args.out #path for the retrieval data to be written to

    decomp_mol = []
    #extract the decomposition data
    for index,line in enumerate(decompo_set_file):
        decomp_data = line.strip("\n").split("|")
        decomp_mol.append(decomp_data)
    decompo_set_file.close()

    #perform building block retrieval
    RB = retrieve_building_blocks(mol_set,output_path,mol_limit)
    RB.retrieve_bb_pairs(decomp_mol)

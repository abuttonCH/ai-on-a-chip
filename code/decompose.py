import numpy as np
from rdkit.Chem import MolToSmiles as MS
from reaction_library import reaction_library
from rdkit.Chem import rdChemReactions as R
import time
import os
import sys
import argparse

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

'''
This code performed the decomposition reactions to produce the corresponding
reactant molecules.
'''

class retro_synthesis:
    def __init__(self,input_filename,output_filepath,limit):
        #loads the input file
        self.input_file = open(input_filename,"r")
        #stores input molecule smiles from the input file
        self.smiles_list = []
        for lines in self.input_file:
            self.smiles_list.append(lines.strip("\n"))
        #writes the output file
        self.output_file = open(output_filepath,"w+")
        self.output_file.write("ID|product molecule|decomposition reaction|decomposed reactants\n")
        #limit the number of molecules that are read in
        self.limit = limit

    def decompose(self,mol_id,smiles,reaction_func):
        #perform the decomposition reaction
        reactant_list = reaction_func.run_reaction(smiles)
        output_list = []
        for reactant_mol in reactant_list:
            reactant_pair = []
            for mol in reactant_mol:
                #converts reactant rdkit mol objects to smiles
                reactant_pair.append(MS(mol))
            #makes sure we don't log duplicates of the decomposition reaction
            if reactant_pair not in output_list:
                output_list.append(str(mol_id)+"|"+smiles+"|"+reaction_func.reaction_name+"|"+".".join(reactant_pair))
        return output_list

    def run_through_a_list(self,smiles_list,smarts_list):
        output = []
        for rxn_idx,smarts_rxn in enumerate(smarts_list):
            for product_index,mol in enumerate(smiles_list):
                #limit the number of molecules that are read in
                if self.limit != None:
                    if product_index > int(self.limit):
                        break
                output = output + self.decompose(product_index,mol,smarts_rxn)
        return output

    def full_decomposition(self,smiles_list,smarts_list):
        decomp_output = []
        reactant_list = self.run_through_a_list(self.smiles_list,smarts_list)
        decomp_output = decomp_output + reactant_list
        while len(reactant_list) != 0:
            reactant_list = self.run_through_a_list(reactant_list,smarts_list)
            decomp_output = decomp_output + reactant_list

        string_output = ".".join(list(np.unique(decomp_output)))
        print("full decomp\n",string_output)

    def single_step_decomposition(self,smarts_list):
        decomp_output = self.run_through_a_list(self.smiles_list,smarts_list)
        #decomp_output = decomp_output + reactant_list
        for decomp_pair in np.unique(decomp_output): #prevents storing duplicates
            self.output_file.write(decomp_pair+"\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mol')
    parser.add_argument('--reaction')
    parser.add_argument('--out')
    parser.add_argument('--limit')
    args = parser.parse_args()

    #inputs
    input_file_path = args.mol #molecules to be decomposed
    reaction_file_path = args.reaction #decomposition reactions
    mol_limit = args.limit #limit on the number of molecules
    output_file_path = args.out #path for the decomposition data to be written to

    #loads the reaction file
    reaction_file = open(reaction_file_path,"r")
    reaction_header = reaction_file.readline()

    #stores the reaction info
    reaction_list = []
    #reads in and separates the reaction data
    for reaction_line in reaction_file:
        Name, SMARTS, ring_change_count = reaction_line.strip("\n").split("|")
        try:
            reaction_list.append(reaction_library(Name,SMARTS,ring_change_count))
        except:
            print(Name,SMARTS)
    reaction_file.close()

    RS = retro_synthesis(input_file_path,output_file_path,mol_limit)
    RS.single_step_decomposition(reaction_list)

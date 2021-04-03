from rdkit.Chem import rdChemReactions
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from rdkit.Chem.rdMolDescriptors import *

class reaction:
    def __init__(self,reaction_name,reaction,reaction_label):
        self.reaction_name = reaction_name
        self.reaction = AllChem.ReactionFromSmarts(reaction)
        rdChemReactions.ChemicalReaction.Initialize(self.reaction)
        self.reaction_label = reaction_label

    def Check_number_reactants(self):
        return self.reaction.GetNumReactantTemplates()

    def IsMoleculeReactant(self,input_mol):
        Is_building_block_a_Reactant = rdChemReactions.ChemicalReaction.IsMoleculeReactant(self.reaction, input_mol)
        condition = None
        if input_mol == condition:
            Is_building_block_a_Reactant = False
        return Is_building_block_a_Reactant

    def Reactant_position(self,mol):
        Reactants_generalized = self.reaction.GetReactants()
        Reactant_positions = []
        for position, position_check in enumerate([Chem.Mol.HasSubstructMatch(mol, reactant) for reactant in Reactants_generalized]):
            if position_check == True:
                Reactant_positions.append(position)
        return Reactant_positions

    def Reaction_components(self):
        Reactants_generalized = self.reaction.GetReactants()
        Reactants_smarts = [Chem.MolToSmiles(mol).replace('*','C') for mol in Reactants_generalized]
        Reactants_components = [Chem.MolFromSmiles(mol) for mol in Reactants_smarts]

        Product_generalized = self.reaction.GetProducts()
        Product_smarts = Chem.MolToSmiles(Product_generalized[0]).replace('*','C')
        Product_components = Chem.MolFromSmiles(Product_smarts)

        return Reactants_components+[Product_components]

    def PerformReaction(self,input_mol):
        product = []
        #apply the inputs in all possible reactant positions
        for shift in range(0,len(input_mol)):
            product = np.append(product,np.concatenate([self.reaction.RunReactants(list(np.roll(input_mol,shift)))]))

        product = self.Product_condition(input_mol,product)
        if len(product) != 0:
            #removes duplicates
            product_smiles = [Chem.MolToSmiles(mol) for mol in product]
            product_smiles = list(set(product_smiles))
            product = [Chem.MolFromSmiles(mol) for mol in product_smiles]

        return product

    def Reactant_condition(self,input_mol):
        if self.reaction_name == "Lactonization":
            for mol in input_mol:
                # Product condition
                pass
        return input_mol

    def Product_condition(self,input_mol,product):
        if self.reaction_name == "Lactonization":
            product_list = []
            for mol in product:
                # Product condition
                try:
                    ring_count = CalcNumRings(input_mol[0])
                    Chem.SanitizeMol(mol[0])
                    ring_change = CalcNumRings(mol[0]) - ring_count
                    if ring_change == 1:
                        # performs the reaction
                        product_list.append(mol)
                except:
                    continue
        elif self.reaction_name == "Ester_formation":
            product_list = []
            for mol in product:
                # Product condition
                try:
                    #checks the net change in the number of carxoblyic acids
                    Chem.SanitizeMol(mol[0])
                    Acid_anhydride_count_input_mol = len(input_mol[0].GetSubstructMatches(Chem.MolFromSmiles("O=COC=O")))
                    Acid_anhydride_count_mol = len(mol[0].GetSubstructMatches(Chem.MolFromSmiles("O=COC=O")))
                    Acid_anhydride_change = abs(Acid_anhydride_count_input_mol - Acid_anhydride_count_mol)
                    if Acid_anhydride_change == 0:
                        # performs the reaction
                        product_list.append(mol)
                except:
                    continue
        elif self.reaction_name == ("FGI_bromination" or "FGI_clorination"):
            product_list = []
            for mol in product:
                # Product condition
                try:
                    Chem.SanitizeMol(mol[0])
                    Carboxlyic_count_input_mol = len(input_mol[0].GetSubstructMatches(Chem.MolFromSmiles("C(=O)O")))
                    Carboxlyic_count_mol = len(mol[0].GetSubstructMatches(Chem.MolFromSmiles("C(=O)O")))
                    Carboxlyic_change = Carboxlyic_count_input_mol - Carboxlyic_count_mol
                    if Carboxlyic_change == 0:
                        # performs the reaction
                        product_list.append(mol)
                except:
                    continue
        #if there is no condition for the products
        else:
            product_list = product

        return product_list

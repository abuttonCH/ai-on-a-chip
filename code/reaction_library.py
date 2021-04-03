import numpy as np
from rdkit.Chem import MolToSmiles as MS
from rdkit.Chem import MolFromSmiles as SM
from rdkit.Chem import rdChemReactions as R
from rdkit.Chem.rdMolDescriptors import CalcNumRings
from rdkit import Chem

class reaction_library:
    def __init__(self,reaction_name,smarts,ring_change_count):
        self.reaction_name = reaction_name
        self.smarts = smarts
        self.rxn = R.ReactionFromSmarts(smarts)
        self.ring_change_count = int(ring_change_count)

    def set_reaction(self,product_smiles):
        try:
            product_mol = SM(product_smiles)
            prod_num_rings = CalcNumRings(product_mol)
        except:
            print("error in the product",product_smiles)
            return []

        try:
            reactant_list = self.rxn.RunReactants([product_mol])
        except:
            print("Reaction failed")
            print(self.reaction_name,self.smarts,product_smiles)
            exit()

        approved_reactants = []
        for reactant_mol in reactant_list:
            #condition 1 - conserved ring count
            try:
                [Chem.SanitizeMol(r) for r in reactant_mol]
                if np.sum([CalcNumRings(r) for r in reactant_mol]) - prod_num_rings == self.ring_change_count:
                    approved_reactants.append(reactant_mol)
            except:
                print("could not sanitize ",product_smiles,".".join([MS(r) for r in reactant_mol]))
        return approved_reactants

    def run_reaction(self,product_smiles):
        return self.set_reaction(product_smiles)

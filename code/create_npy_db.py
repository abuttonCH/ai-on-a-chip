import numpy as np
import argparse

def create_db(mol_input,mol_output):
    mol_db = open(mol_input,"r")
    header = mol_db.readline()

    mol_db_array = []
    for line in mol_db:
        mol_data = line.strip("\n").split(",")
        mol_db_array.append(mol_data)

    mol_db_array = np.asarray(mol_db_array,dtype=bytes)
    np.save(mol_output,mol_db_array)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input')
    parser.add_argument('--output')
    args = parser.parse_args()

    create_db(args.input,args.output)

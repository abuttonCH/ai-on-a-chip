#create the mol_db numpy array object
echo "creating mol_db"
python code/create_npy_db.py --input data/mol_db_data.csv --output data/LSTM_FLOW-MOL_DB_DATA.npy
echo "done"
#convert the molecules into their corresponding reactants
echo "converting molecules to reactants"
python code/decompose.py --mol data/data_val.txt --reaction data/decomposition_reactions.txt --out output/test_decomp.txt --limit 100
echo "done"
#compare the reactants to the mol_db molecules
echo "retrieving building block molecules"
python code/retrieve_bb.py --decomp output/test_decomp.txt --mol_db data/LSTM_FLOW-MOL_DB_DATA.npy --out output/test_retrieve.txt
echo "done"

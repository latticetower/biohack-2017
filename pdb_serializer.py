
import pickle 
import argparse

PDB_PATH = "pdb"

import os, prody

if not os.path.exists(PDB_PATH):
    os.mkdir(PDB_PATH)

prody.proteins.localpdb.pathPDBFolder(PDB_PATH)

def downloadPDBStructures(pdbIds):
    """
	download all structures if they are not present and were not previously processed 
	returns list of structures which were actually downloaded
	"""
    downloadedPDBs = zip(pdbIds, prody.proteins.localpdb.fetchPDB(pdbIds)) 
    downloadedPDBs = map(lambda x: x[0], filter(lambda x: x[1] is not None, downloadedPDBs)) 
    return downloadedPDBs
	
if __name__ == "__main__":
	gene_info_file = "gene_name_to_gene_ids_pdb_ids.pickle"
	if not os.path.exists(gene_info_file):
		print("File %s couldn't be found" % gene_info_file)
		exit(1)
	data = pickle.load(open(gene_info_file,"rb"))
	# data.values() - contains (geneIds, pdbIds)
	set_pdb = set()
	for (_, pdbIds) in data.values():
		set_pdb |= set(downloadPDBStructures(pdbIds))
	# just in case we want to save the results

	with open('set_pdb.pickle', 'wb') as f:
		pickle.dump(set_pdb, f)


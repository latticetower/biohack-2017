from Bio.PDB import *
import numpy as np

def get_cys_atoms(model, chain):
    cys_atoms = []
    for residue in model[chain]:
        if residue.resname == 'CYS':
            for atom in residue:
                if atom.get_name() == 'SG':
                    cys_atoms.append(atom)
                    #print(residue.resname)
    return cys_atoms
            
def get_calpha_atoms(model, chain):
    return [atom for atom in model[chain].get_atoms()
                     if atom.get_name() == 'CA']


def find_min_dist(cas,sgs):
    result = []
    min = 100
    for sg in sgs:
        for ca in cas:
            dist = (sg.get_vector()-ca.get_vector()).norm()
            if dist < 10:
                if dist < min:
                    min = dist
                    result = [(dist, ca.get_parent().id[1])]
    return result

def get_coords(atoms):
    return np.array([atom.get_coord() for atom in atoms])

def find_sites(pdb_data, model):
    result = []

    pdb_id = pdb_data[0]
    chain_oncogene = pdb_data[1]
    chain_petide = pdb_data[2]

    calpha_atoms = get_calpha_atoms(model, chain_petide)
    cys_atoms = get_cys_atoms(model, chain_oncogene)

    distances = find_min_dist(calpha_atoms, cys_atoms)

    for tup in distances:
        result.append((pdb_id, chain_oncogene, chain_petide, tup[1], tup[0]))

    return result


def find_all_sites(pdb_info):
    result = []
    for tup in pdb_info:
        parser = PDBParser()
        structure = parser.get_structure('test', tup[0])
        model = list(structure.get_models())[0]
        for ch in model:
            if tup[1] != ch.id:
                result += find_sites((tup[0], tup[1], ch.id), model)
    return result


# pdb_info = [('3ny5.pdb', ['A','B','C','D'], ['A','B','C','D']), ('5fcg.pdb', ['A'], ['C'])]
pdb_info = [('5fcg.pdb', 'A')]
print(find_all_sites(pdb_info))
#('3ny5.pdb', 'A'), ('3ny5.pdb', 'B'), ('3ny5.pdb', 'C'), ('3ny5.pdb', 'D'),
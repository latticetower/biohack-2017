from bioservices.kegg import KEGG
k = KEGG()

import pickle
import argparse

ORGANISM = "hsa"
GENES = ["p53"] # sample gene
PDB_PATH = "pdb"

import os, prody, pystache

if not os.path.exists(PDB_PATH):
    os.mkdir(PDB_PATH)
    # TODO: for now I haven't checked if pathPDBFolder creates this folder - 
    # if it is created, this check should be removed.
prody.proteins.localpdb.pathPDBFolder(PDB_PATH)


def getKeggInfoForGene(organism, gene_name):
        geneIdsInPathways = [ 
            x.split()[0]
            for x in keggParser.find(organism, gene_name).split("\n") 
            if len(x.split()) > 0 ]
        geneIdByPDB = dict()
        # for now store sequence information for gene in dictionary also
        uniprotIdsByGene = dict()
        for geneId in geneIdsInPathways:
            geneInfo = keggParser.parse(keggParser.get(geneId))
            if 'DBLINKS' in geneInfo:
                if 'UniProt' in geneInfo['DBLINKS']:
                    uniprotIdsByGene[geneId] = set(geneInfo['DBLINKS']['UniProt'].split())
                # this might be empty, should show something in else clause
            if 'STRUCTURE' in geneInfo:
                for pdb_name in geneInfo['STRUCTURE']['PDB'].split():
                    geneIdByPDB[pdb_name] = geneId
        pdbIds = list(geneIdByPDB)
        return pdbIds, geneIdByPDB, uniprotIdsByGene


def getGeneAssociatedPDBStructures(organism, gene_name):
    (pdbIds, geneIdByPDB, uniprotIdsByGene) = getKeggInfoForGene(organism, gene_name)

    # next - download all structures if they are not present and were not previously processed 

    # all other checks are skipped now, will appear here later
    downloadedPDBs = zip(pdbIds, prody.proteins.localpdb.fetchPDB(pdbIds)) # use filtered structures list
    # TODO: some of PDB identifiers correspond to structures which has only .cif files. 
    # in case if they are present we should check fetched PDB files and to store PDBIDs 
    # which were not downloaded somewhere.
    downloadedPDBs = map(lambda x: x[0], filter(lambda x: x[1] is not None, downloadedPDBs)) # this filters zero values
    return downloadedPDBs



def get_reference_chain(pdbid, geneIds):
    """
    returns chain 
    """
    uniprotIdDict = geneIds
    atoms_header = prody.parsePDBHeader(pdbid)
    for polymer in atoms_header['polymers']:
        uniprotRef = filter(lambda k: keggParser.database == 'UniProt', 
            atoms_header[polymer.chid].dbrefs)
        if len(uniprotRef) < 1:
            continue
        if len(set(map(lambda x: x [0].accession, uniprotRef)) & uniprotIdDict) > 0:
            return polymer.chid
    return None


def iterate_over_pairs(pdbid, ref_chid):
    """
    iter<ref_chain, other_chain>
    """
    header = prody.parsePDBHeader(pdbid)
    for p in header['polymers']:
        if p.chid == ref_chid:
            continue
        yield (ref_chid, p.chid)


  

def getPairInformation(pdbid, reference_chain, pair_chain, cutoff=5, covalent_bond_cutoff=5):
    """
    1. reads pdb id from file
    2. selects atoms from pair of chains within cutoff
    # draws selection within interface (simplest possible view)
    reference == oncogene,
    pair == peptide
    """
    atoms = prody.parsePDB(pdbid) # TODO: turn off debug
    reference_atoms = atoms.select("chain %s and not water" % reference_chain)
    pair_atoms = atoms.select("chain %s and not water" % pair_chain)
    # next try to select everything
    ref_contacts = prody.measure.contacts.Contacts(reference_atoms)
    ref_selection = ref_contacts.select(cutoff, pair_atoms) # we need these atoms
    print(ref_selection)
    pair_contacts = prody.measure.contacts.Contacts(pair_atoms)
    pair_selection = pair_contacts.select(cutoff, reference_atoms) # and these
    sulfur_pairs = []
    ## 1. select Cys atoms on oncogene
    for (r, ch2, distance) in prody.measure.contacts.findNeighbors(reference_atoms, covalent_bond_cutoff, pair_atoms):
        if r.getResname() in 'CYS': # and r.getElement() in ['S'] :
            sulfur_pairs.append((r.getSerial(), ch2.getSerial()))
    # filtering: if there is no Cys, return nothing
    if len(sulfur_pairs) < 1:
        return None
    return (pdbid, reference_chain, set(ref_selection.getResnums()), 
        pair_chain, set(pair_selection.getResnums()), sulfur_pairs)    
    #prody.proteins.functions.showProtein(reference_atoms, pair_atoms);


def check_condition(pdbid, ref_chid, other_chid):
    """
    main condition for selection pair of chains based on interactions
    return boolean value and interface information (to draw later)
    """
    #TODO: obsolete now. use getPairInformation instead
    (pdbid, ref, res_ref_no, pair_chain, ref, sulfur_pairs) = getPairInformation(pdbid, ref_chid, other_chid) 
    if len(ref_chid) < 0 and len(ref_ch2) < 0:
        return False, (p, r, no, pc, sr, rch)
    if len(sulfur_pairs) < 0: return None
    #else return 
    pass 


def get_candidates_for(pdbid, ref_chid):
    """
    Selects good pairs for ref_chid based on set of parameters
    """
    for (ref, other) in iterate_over_pairs(pdbid, ref_chid):
        condition, interface_info = check_condition(pdbid, ref, other)
        if condition:
            yield other, interface_info


def iterate_over_objects():
    # 1. 
    gene_info_file = "gene_name_to_gene_ids_pdb_ids.pickle"
    if not os.path.exists(gene_info_file):
        print("File %s couldn't be found" % gene_info_file)
        return exit(1)
    data = pickle.load(open(gene_info_file,"rb"))
    set_pdb =  pickle.load(open('set_pdb.pickle', 'rb'))
    
    for gene_name in data:
        (geneIds, pdbIds) = data[gene_name] 
        downloadedPDBs = pdbIds & set_pdb
        for pdbid in downloadedPDBs:
            chid = get_reference_chain(pdbid, geneIds)
            if chid is None:
                continue # skip if there is no UniProt
            for (other, interface_info) in get_candidates_for(pdbid, chid):
                yield other
            #print(chid)
    

def renderTemplate(template, info):
    """
    saves pymol script to pictures.
    if run from current folder, this scripts saves pdb id to file.
    """
    renderer = pystache.Renderer()
    lines = []
    with open(template, 'r') as k:
        for i  in k.readlines():
            (pdbid, reference, ref_selection, pair_chain, pair_selection, sulfur_pairs) = info
            s = renderer.render(i, {
                'pdbid': pdbid,
                'onco_chain': reference,
                'peptide_chain': pair_chain,
                'onco_resnum': "+".join(map(str, ref_selection)), 
                'peptide': '+'.join(map(str, pair_selection)),
                'pept_cys': '+'.join(map(lambda x: str(x[1]), sulfur_pairs))
                })
            lines.append(s)
    if not os.path.exists("pictures"):
        os.mkdir("pictures")
    image = os.path.join("pictures"+reference+"_"+pair_chain+".pml")
    with open(image, 'w') as f:
        for line in lines:
            f.write(line)
    cmd([" ".join(["pymol -c", "pictures"+reference+"_"+pair_chain+".pml"])])
    pass


from subprocess import call

if __name__== "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pdbid", help="PDBID", type=str)
    parser.add_argument("oncogene_chain", help="oncogene chain", type=str)
    #parser.add_argument("peptide_chain", help="peptide chain", type=str)
    args = parser.parse_args()
    for other in iterate_over_pairs(pdbid, oncogene_chain):
        info = getPairInformation(pdbid, oncogene_chain, other)
        #info = getPairInformation("2OSL", "L", "H")
        renderTemplate("structure_view.pml.mustache", info)

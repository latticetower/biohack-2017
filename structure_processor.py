"""
If run:

```
python3 structure_processor.py "" "" --filter_genes "TP53"
```
saves to different files data for gene named "TP53" (this parameter can be comma-separated list of gene names).
Saves to pictures/ fragments if they are found.
Otherwise saves to different files.
after processing saves to "processed_genes.log" gene names from parameter list. 
To rerun with the same gene list, remove lines corresponding to names from this file or remove the whole file - 
currently it is used to skip gene names which were already processed.

Another call option might be incorrect now.
"""
from bioservices.kegg import KEGG
keggParser = KEGG()

import pickle
import argparse

ORGANISM = "hsa"
GENES = ["p53"] # sample gene
PDB_PATH = "pdb"

import os, prody, pystache, logging

if not os.path.exists(PDB_PATH):
    os.mkdir(PDB_PATH)
    # TODO: for now I haven't checked if pathPDBFolder creates this folder - 
    # if it is created, this check should be removed.
prody.proteins.localpdb.pathPDBFolder(PDB_PATH)

keggCache = dict()
cache_file_name = "kegg_info_cache.pickle"
if os.path.exists(cache_file_name):
    keggCache = pickle.load(open(data_file, "rb"))

def getKeggInfoForGene(organism, gene_name):
        if gene_name  in keggCache:
            return keggCache[gene_name]
        geneIdsInPathways = [ 
            x.split()[0]
            for x in keggParser.find(organism, gene_name).split("\n") 
            if len(x.split()) > 0 ]
        geneIdByPDB = dict()
        # for now store sequence information for gene in dictionary also
        uniprotIdsByGene = dict()
        pdbIdsByGene = dict()
        aaSequencesByGene = dict()
        for geneId in geneIdsInPathways:
            geneInfo = keggParser.parse(keggParser.get(geneId))
            if 'DBLINKS' in geneInfo:
                if 'UniProt' in geneInfo['DBLINKS']:
                    uniprotIdsByGene[geneId] = set(geneInfo['DBLINKS']['UniProt'].split())
                # this might be empty, should show something in else clause
            if 'STRUCTURE' in geneInfo:
                for pdb_name in geneInfo['STRUCTURE']['PDB'].split():
                    geneIdByPDB[pdb_name] = geneId
                pdbIdsByGene[geneId] = geneInfo['STRUCTURE']['PDB'].split()
            if 'AASEQ' in geneInfo:
                aaSequencesByGene[geneId] = geneInfo['AASEQ'].replace(' ', '')
        pdbIds = list(geneIdByPDB)
        keggCache[gene_name] = (pdbIds, geneIdByPDB, uniprotIdsByGene, aaSequencesByGene, pdbIdsByGene)
        with open(cache_file_name, "wb") as f:
            pickle.dump(keggCache, f)  
        return pdbIds, geneIdByPDB, uniprotIdsByGene, aaSequencesByGene, pdbIdsByGene

prody.LOGGER.verbosity = "none"

def getGeneAssociatedPDBStructures(organism, gene_name):
    (pdbIds, geneIdByPDB, uniprotIdsByGene) = getKeggInfoForGene(organism, gene_name)

    # next - download all structures if they are not present and were not previously processed 
    pdbs = prody.proteins.localpdb.fetchPDB(pdbIds)
    if isinstance(pdbs, str):
        pdbs = [pdbs]
    # all other checks are skipped now, will appear here later
    downloadedPDBs = zip(pdbIds, pdbs) # use filtered structures list
    # TODO: some of PDB identifiers correspond to structures which has only .cif files. 
    # in case if they are present we should check fetched PDB files and to store PDBIDs 
    # which were not downloaded somewhere.
    downloadedPDBs = list(map(lambda x: x[0], filter(lambda x: x[1] is not None, downloadedPDBs))) # this filters zero values
    return downloadedPDBs


def downloadPdbList(pdbIds):
    if isinstance(pdbIds, set):
        return downloadPdbList(list(pdbIds))
    # next - download all structures if they are not present and were not previously processed 
    pdbs = prody.proteins.localpdb.fetchPDB(pdbIds)
    if isinstance(pdbs, str):
        pdbs = [pdbs]
    if pdbs is None:
        pdbs = []
    # all other checks are skipped now, will appear here later
    downloadedPDBs = zip(pdbIds, pdbs) # use filtered structures list
    # TODO: some of PDB identifiers correspond to structures which has only .cif files. 
    # in case if they are present we should check fetched PDB files and to store PDBIDs 
    # which were not downloaded somewhere.
    downloadedPDBs = list(map(lambda x: x[0], filter(lambda x: x[1] is not None, downloadedPDBs))) # this filters zero values
    return downloadedPDBs

def get_reference_chain(pdbid, geneIds):
    """
    returns chain 
    """
    uniprotIdDict = geneIds
    try:
        atoms, atoms_header = prody.parsePDB(pdbid, header=True)
    except:
        with open("problem_with_parsing.txt", 'a') as f:
            f.write(pdbid+"\n")
        return None
    #print(list(atoms_header))
    for polymer in atoms_header['polymers']:
        uniprotRef = [ r.accession
            for r in atoms_header[polymer.chid].dbrefs
            if r.database == 'UniProt'
            ]
        if len(uniprotRef) < 1:
            continue
        # print(uniprotRef, uniprotIdDict)
        if len(set(uniprotRef) & uniprotIdDict) > 0:
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


  

def getPairInformation(pdbid, reference_chain, pair_chain, cutoff=5, covalent_bond_cutoff=2.5):
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


def check_condition(pdbid, ref_chid, other_chid, bonds1=5, bonds2=2.5):
    """
    main condition for selection pair of chains based on interactions
    return boolean value and interface information (to draw later)
    """
    #TODO: obsolete now. use getPairInformation instead
    (pdbid, ref, res_ref_no, pair_chain, pair_res_no, sulfur_pairs) = getPairInformation(
        pdbid, ref_chid, other_chid, bonds1, bonds2) 
    if len(res_ref_no) < 0 and len(pair_res_no) < 0:
        print(111)
        return False, (pdbid, ref, res_ref_no, pair_chain, pair_res_no, sulfur_pairs)
    if len(sulfur_pairs) < 0:
        return False, None
    return True, (pdbid, ref, res_ref_no, pair_chain, pair_res_no, sulfur_pairs)
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
    #TODO: obsolete
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
                pass
                #yield other
            #print(chid)
    

def renderTemplate(template, info):
    """
    saves pymol script to pictures.
    if run from current folder, this scripts saves pdb id to file.
    """
    renderer = pystache.Renderer()
    lines = []
    (pdbid, reference, ref_selection, pair_chain, pair_selection, sulfur_pairs) = info
    with open(template, 'r') as k:
        for i  in k.readlines():
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
    image = os.path.join("pictures", reference+"_"+pair_chain+".pml")
    with open(image, 'w') as f:
        for line in lines:
            f.write(line)
    cmd([" ".join(["pymol -c", "pictures/"+reference+"_"+pair_chain+".pml"])])
    pass


def appendToNotebook(template, info, gene_name):
    """
    saves pymol script to pictures.
    if run from current folder, this scripts saves pdb id to file.
    """
    renderer = pystache.Renderer()
    lines = []
    (pdbid, reference, ref_selection, pair_chain, pair_selection, sulfur_pairs) = info
    with open(template, 'r') as k:
        for i  in k.readlines():
            s = renderer.render(i, {
                'pdbid': pdbid,
                'onco_chain': reference,
                'peptide_chain': pair_chain,
                'onco_resnum': "+".join(map(str, ref_selection)), 
                'peptide': '+'.join(map(str, pair_selection)),
                'pept_cys': '+'.join(map(lambda x: str(x[1]), sulfur_pairs))
                })
            lines.append(s)
    if not os.path.exists("notebooks"):
        os.mkdir("notebooks")
    image = os.path.join("notebooks", gene_name+".ipynb")
   
    with open(image, 'a') as f:
        for line in lines:
            f.write(line)
    #cmd([" ".join(["pymol -c", "pictures"+reference+"_"+pair_chain+".pml"])])
    pass

def doFilter(geneNames, b1=5, b2=2.5):
    genesToReturn = []
    if os.path.exists("processed_genes.log"):
        with open("processed_genes.log") as f:
            processedGenes = set([line.strip() for line in f])
    else:
        processedGenes=set()
    
    gene_info_file = "gene_name_to_gene_ids_pdb_ids_filtered_by_cys.pickle"
    if not os.path.exists(gene_info_file):
        print("File %s couldn't be found" % gene_info_file)
        return
    genes = pickle.load(open(gene_info_file, "rb"))
    set_pdb =  pickle.load(open('set_pdb.pickle', 'rb'))
    
    for gene_name in geneNames:
        if gene_name in processedGenes:
            continue
        #print(gene_name)
        (geneIds, pdbIds) = genes[gene_name]
        downloaded = downloadPdbList(pdbIds)
        #print(len(downloaded))
        with open("unavailable_pdb_list.txt", 'a') as f:
            for pdbid in (set(pdbIds) - set(downloaded)):
                f.write("%s\n"% pdbid)
        
        with open("downloaded_pdbs.txt", 'a') as g:
            g.write(gene_name + " " + " ".join(downloaded)+"\n")
        
        has_results = False
        for pdbid in downloaded:
            chid = get_reference_chain(pdbid, geneIds)
            if chid is None:
                with open("bad_pdb_for_gene_list.txt", 'a') as g:
                    g.write(gene_name + " " + pdbid + "\n")
                continue # skip if there is no UniProt

            no_contacts = True
            for (ref, other) in iterate_over_pairs(pdbid, chid):
                condition, interface_info = check_condition(pdbid, ref, other, b1, b2)
                if condition:
                    no_contacts = False
                    has_results=True
                    with open("ppi_list.txt", 'a') as f:
                        f.write(" ".join([
                        pdbid, ref, other, str(b1), str(b2)
                        ])+"\n")
                    renderTemplate("structure_view.pml.mustache", interface_info)
                    appendToNotebook("notebook.ipynb.mustache", interface_info, gene_name)
            if no_contacts:
                with open("no_contacts_with_gene.txt", 'a') as g:
                    g.write(gene_name + " " + pdbid + " " + chid+ "\n")
        if has_results:
            print("gene", gene_name, "has some Cys related things")
            genesToReturn.append(gene_name)
        else: 
            print("no results for gene", gene_name)           
        # add to already processed list
        #logging.debug(gene_name)
        with open("processed_genes.log", 'a') as f:
            f.write(gene_name+"\n")
    if len(genesToReturn) > 0:
        return genesToReturn
    pass
    
    
from subprocess import call
logging.basicConfig(filename='processed_genes.log', format='%(message)s', level=logging.DEBUG)


def parse():
    """temporary to save kegg values"""
    gene_info_file = "gene_name_to_gene_ids_pdb_ids_filtered_by_cys.pickle"
    if not os.path.exists(gene_info_file):
        print("File %s couldn't be found" % gene_info_file)
        return exit(1)
    genes = pickle.load(open(gene_info_file, "rb"))
    keggIds = list(set([y for x in genes.values() for y in x[0]]))
    i = 0
    while i < len(keggIds):
        res = keggParser.get("+".join(keggIds[i: min(i+20, len(keggIds))]))
        i+=20
        with open("keggoutput.txt", 'a') as f:
            f.write(res+"\n")
    #print(kegg.parse(res))
    pass
    
def process():
    entry = None
    ids = set()
    inDBLINKS = False
    keggToUniprot = dict()
    with open("keggoutput.txt") as f:
        for line in f:
            if line.startswith("ENTRY"):
                entry = "hsa:"+line.split()[1]
                if inDBLINKS and entry is not None:
                    keggToUniprot[entry] = list(ids)
                    entry = None
                continue
            if line.startswith("DBLINKS"):
                if line.strip().split()[1].startswith("UniProt"):
                    ids |= set(line.strip().split()[2:])
                inDBLINKS = True
                continue
            if line.startswith(' '):
                if inDBLINKS:
                    #print(line.strip().split()[0])
                    if line.strip().split()[0].startswith("UniProt"):
                        print(line.strip().split()[1:])
                        ids |= set(line.strip().split()[1:])
                    pass
            else:
                inDBLINKS = False
                if entry is not None:
                    keggToUniprot[entry] = list(ids)
                    entry = None
    with open("keggToUniprot.pickle", "wb") as f:
        pickle.dump(keggToUniprot, f)
    #print( keggToUniprot)
        
            
    
if __name__== "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pdbid", help="PDBID", type=str)
    parser.add_argument("oncogene_chain", help="oncogene chain", type=str)
    parser.add_argument("--filter_genes", default="")
    #parser.add_argument("peptide_chain", help="peptide chain", type=str)
    args = parser.parse_args()
    if len(args.filter_genes)>0:
        doFilter(args.filter_genes.split(","))
        exit()
        
    for other in iterate_over_pairs(args.pdbid, args.oncogene_chain):
        info = getPairInformation(args.pdbid, args.oncogene_chain, other)
        #info = getPairInformation("2OSL", "L", "H")
        renderTemplate("structure_view.pml.mustache", info)

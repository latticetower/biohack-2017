from bioservices import KEGG
import os

def kegg_find(*args):
    if not hasattr(kegg_find,"cache"):
        if os.path.isfile("kegg_find.cache"):
            kegg_find.cache = pickle.load(open("kegg_find.cache","rb"))
        else:
            kegg_find.cache = {}

    if args not in kegg_find.cache or kegg_find.cache[args] is None:
        k = KEGG()
        result = k.find(*args)
        kegg_find.cache[args] = result
        with open("kegg_find.cache~","wb") as f:
            pickle.dump(kegg_find.cache, f)
        os.rename("kegg_find.cache~", "kegg_find.cache")
        return result
    else:
        return kegg_find.cache[args]

def kegg_get(*args):
    if not hasattr(kegg_get,"cache"):
        if os.path.isfile("kegg_get.cache"):
            kegg_get.cache = pickle.load(open("kegg_get.cache","rb"))
        else:
            kegg_get.cache = {}

    if args not in kegg_get.cache or kegg_get.cache[args] is None:
        k = KEGG()
        result = k.get(*args)
        kegg_get.cache[args] = result
        with open("kegg_get.cache~","wb") as f:
            pickle.dump(kegg_get.cache, f)
        os.rename("kegg_get.cache~", "kegg_get.cache")
        return result
    else:
        return kegg_get.cache[args]


def get_pdb_id_by_name_gene(gene_name):
    k = KEGG()
    gene_ids = []

    # по названию гена получаем id
    gen = kegg_find("hsa", gene_name)
    if gen in [400, 404]:
        return [],[]
    for line in gen.split("\n"):
        if len(line)>0:
            gene_ids.append(line.split("\t")[0])

    # по каждому полученному id гена получаем PDB_ID
    pdb_ids = []
    if len(gene_ids)>100:
        return [],[]
    for gene_id in gene_ids:

        e = kegg_get(gene_id)
        if e in [400, 404]:
            continue
        d = k.parse(e)

        if "STRUCTURE" in d:
            pdb_ids += d["STRUCTURE"]["PDB"].split()
    return gene_ids, pdb_ids

gene_name_to_gene_ids_pdb_ids = {}

if __name__ == "__main__":
    import pickle
    genes = pickle.load(open("oncogen_genes.pickle","rb"))
    print(genes)
    for i, gene in enumerate(genes):
        print(i, gene)
        gene_ids, pdb_ids = get_pdb_id_by_name_gene(gene)
        print("%2d / %2d "%(i, len(genes)), gene, "{%d}{%d}"%(len(gene_ids),len(pdb_ids)))
        gene_name_to_gene_ids_pdb_ids[gene] = get_pdb_id_by_name_gene(gene)
    with open("gene_name_to_gene_ids_pdb_ids.pickle","wb") as f:
        pickle.dump(gene_name_to_gene_ids_pdb_ids,f)

from bioservices import KEGG


def list_of_genes(genes):
    all = []
    for gene in genes:
        all.append(get_pdb_id_by_name_gene(gene))
    return all


def get_pdb_id_by_name_gene(name_gene):
    k = KEGG()
    id_gene = []

    # по названию гена получаем id
    gen = k.find("hsa", name_gene)
    for i in gen.split():
        if "hsa" in i:
            id_gene.append(i.split(':')[1])

    # по каждому полученному id гена получаем PDB_ID
    pdb_id = []
    for key, value in enumerate(id_gene):
        e = k.get('hsa:{}'.format(id_gene[key]))
        d = k.parse(e)
        if "STRUCTURE" in d:
            pdb_id += d["STRUCTURE"]["PDB"].split()

    return (name_gene, id_gene, pdb_id)

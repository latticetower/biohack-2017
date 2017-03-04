import csv

# returns list of gene names from uniprot
# http://www.uniprot.org/uniprot/?query=keyword:KW-0656
def get_uniprot_oncogenes():
    with open("data/uniprot.tsv") as f:
        rows = list(csv.reader(f, delimiter="\t"))
        column_names = { name: i for i,name in enumerate(rows[0]) }
        rows = rows[1:]
        at = lambda row, column_name : row[column_names[column_name]]
        gene_names = []
        for row in rows:
            gene_names += at(row,"Gene names").split(" ")
        return gene_names


# returns gene_names from Network of cancer genes database
# http://ncg.kcl.ac.uk/download.php
def get_ncg_oncogenes():
    with open("data/ncg.tsv") as f:
        rows = list(csv.reader(f, delimiter="\t"))
        column_names = { name: i for i,name in enumerate(rows[0]) }
        rows = rows[1:]
        at = lambda row, column_name : row[column_names[column_name]]
        gene_names = []
        print(column_names)
        for row in rows:
            if (
                at(row, "potential_false_positive") == "FALSE" and
                at(row,"cancer") == "TRUE"
            ):
                gene_names += [ at(row,"symbol") ]
        return gene_names

if __name__ == "__main__":
    import pickle

    l1 = get_uniprot_oncogenes()
    l2 = get_ncg_oncogenes()
    with open("oncogen_genes.pickle","wb") as f:
        pickle.dump(list(set(l1+l2)), f)



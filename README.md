# biohack-2017
Our team's brand new BioHack (as in biohack.ru) project repository

Steps:

Run in following order:

```
python make_gene_list.py

python kegg.py
```

First command gets a list of gene names from uniprot and entrez databases entries and serializes it.
Second command gets serialized result of first command run and for each gene name 
collects KEGG gene identifiers list and PDB structures list.

`pythin pdb_serializer.py`
- this loads all pdb files and saves all that were correctly downloaded to set.
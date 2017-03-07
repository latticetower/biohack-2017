# biohack-2017
Our team's brand new BioHack (as in biohack.ru) project repository

Our team worked on small research project for Laboratory of Biomolecular NMR (St.Petersburg State University, St. Petersburg, Russia).

The goal of our project was to find a subset of oncogenic target proteins with cysteine located on protein-peptide interaction site. 

At first, we were to detect oncogenic proteins and their interactions from KEGG database.
We retrieved associated Protein Data Bank identifiers for future analysis.

After then, we detected protein-protein interaction interface between oncogene-associated chain and every other chain in PDB file. We also checked for cysteine to be present at the protein surface of the interaction site. If we could find such pair of chains, we saved information about it for future visualization and analysis.

During hackathon we could analyze ~2000 PDB structures out of 13000, and could find ~100 structures with cysteine located at protein-protein interaction interface.


Steps:

Run in following order:

```
python make_gene_list.py

python kegg.py
```

First command gets a list of gene names from uniprot and entrez databases entries and serializes it.
Second command gets serialized result of first command run and for each gene name 
collects KEGG gene identifiers list and PDB structures list.

```
pythin pdb_serializer.py
```
- this loads all pdb files and saves all that were correctly downloaded to set.
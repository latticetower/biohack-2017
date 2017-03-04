import csv
import requests
import io
#pdb_id = '3ny5' #'5fcg #3ny5
#gene = 'BRAF'

def chain_gene(padb_id, gene):
	pattern1 = 'http://www.rcsb.org/pdb/rest/customReport.csv?pdbids='
	pattern2 = '&customReportColumns=geneName&format=csv&service=wsfile'

	r = requests.get(pattern1+pdb_id+pattern2)

	page = list(csv.reader(io.StringIO(r.text)))
	#page = [['structureId', 'chainId', 'geneName'], ['3NY5', 'A', 'BRAF#BRAF1#RAFB1'], ['3NY5', 'B', 'BRAF#BRAF1#RAFB1'], ['3NY5', 'C', 'BRAF#BRAF1#RAFB1'], ['3NY5', 'D', 'BRAF#BRAF1#RAFB1']]
	number_of_chains = len(page) - 1
	chains_gene = [] 
	for line in page:
		line = str(line).split(', ')
		if (line[2].find(gene)!=-1):
			chains_gene += line[1][1]
	dict = {}
	dict['chains_gene'] = chains_gene
	dict['number_of_chains'] = number_of_chains
	return(dict)
print(chain_gene(pdb_id, gene))
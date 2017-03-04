import csv
import requests
import io

def gene_chains(pdb_id, gene):
    pattern1 = 'http://www.rcsb.org/pdb/rest/customReport.csv?pdbids='
    pattern2 = '&customReportColumns=geneName&format=csv&service=wsfile'

    r = requests.get(pattern1+pdb_id+pattern2)

    page = list(csv.reader(io.StringIO(r.text)))
    number_of_chains = len(page) - 1
    chains = []
    print(page)
    for line in page:
        if gene in line[2].split("#"):
            chains += line[1]
    return chains



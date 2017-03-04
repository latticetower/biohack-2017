import subprocess


def calculate_area(pdbid, chain_onco, chain_peptide, cystein_resid):
    subprocess.Popen(['python', 'area.py', pdbid, chain_onco, chain_peptide, cystein_resid])

if __name__ == '__main__':
    import sys
    pdbid, chain_onco, chain_peptide, cystein_resid = sys.argv[1:]
    calculate_area(pdbid, chain_onco, chain_peptide, cystein_resid)

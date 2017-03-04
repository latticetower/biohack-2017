def area(pdbid, chain_onco, chain_peptide, cystein_resid):
    import __main__
    __main__.pymol_argv = ['pymol','-qc']
    import pymol
    from pymol import cmd, stored

    pymol.finish_launching()

    cmd.set('dot_solvent', 3)
    cmd.set('dot_density', 3)

    cmd.load(pdbid)

    cmd.remove('! chain {}+{}'.format(chain_onco, chain_peptide))
    area1 = cmd.get_area('chain {} & resi {} & name SG'.format(chain_onco, cystein_resid))

    cmd.remove('! chain {}'.format(chain_onco))
    area2 = cmd.get_area('chain {} & resi {} & name SG'.format(chain_onco, cystein_resid))

    return area1, area2

if __name__ == '__main__':
    import sys
    pdbid, chain_onco, chain_peptide, cystein_resid = sys.argv[1:]
    first, second = area(pdbid, chain_onco, chain_peptide, cystein_resid)

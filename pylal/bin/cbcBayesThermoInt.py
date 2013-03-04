#!/usr/bin/env python

from argparse import ArgumentParser
import matplotlib.pyplot as pp
import numpy as np

def extract_temp(filename):
    with open(filename, 'r') as inp:
        for line in inp:
            line=line.lstrip().split()
            try:
                i=line.index('Tchain')
                #####################################################
                # WARNING: hardcoded off-by-one adjustment to account
                # for 'null likelihood' header name that splits into
                # two list elements
                #####################################################
                return float(inp.next().split()[i-1])
            except ValueError:
                pass

        raise ValueError('extract_temp: did not find header line with \'Tchain\'')

def get_mean_logl(filename):
    with open(filename, 'r') as inp:
        for line in inp:
            line=line.lstrip().split()
            try:
                line.index('cycle')
                break
            except ValueError:
                pass

        data=np.loadtxt(inp, dtype=[(n, np.float) for n in line])
        N=data.shape[0]

        return np.mean(data['logl'][N/2:])

def thermo_integrands(logls, betas):
    inds=np.argsort(betas)[::-1]

    sorted_logls=logls[inds]
    sorted_betas=betas[inds]

    dbetas=np.diff(np.concatenate((sorted_betas, [0.0])))

    return sorted_betas, -sorted_logls*dbetas

if __name__=='__main__':
    parser=ArgumentParser(description='Thermodynamic integration on PTMCMC samples.')
    
    parser.add_argument('files', metavar='OUTPUT_FILE', nargs='+', help='PTMCMC output files')
    parser.add_argument('--plotfile', metavar='FILE', default='evidence-integrand.pdf', help='plot output file')
    parser.add_argument('--evfile', metavar='FILE', default='evidence.dat', help='evidence output file')

    args=parser.parse_args()

    betas = np.array([1.0/extract_temp(f) for f in args.files])
    logls = np.array([get_mean_logl(f) for f in args.files])

    betas, integrands = thermo_integrands(logls, betas)

    evidence = np.sum(integrands)

    pp.plot(betas, integrands/evidence)
    pp.xscale('log')
    pp.xlabel(r'$\beta$')
    pp.ylabel(r'$\left \langle \log \mathcal{L} \right\rangle d\beta$')
    pp.savefig(args.plotfile)

    with open(args.evfile, 'w') as out:
        out.write('# log(Ev)\n')
        out.write(str(evidence))
    
    
    

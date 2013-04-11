#!/usr/bin/env python

from argparse import ArgumentParser
import matplotlib.pyplot as pp
import numpy as np

def extract_temp(filename):
    """Extracts the PTMCMC temperature from the header lines of the
    given file."""
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
    """Returns the mean value of log(L) from the given filename,
    excluding the first 50% of samples as burn-in."""
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
    """Returns arrays of betas, <log(L)>*d(beta), <log(L)>*d(beta2)
    for the evidence calculation.  The log-evidence is given by
    sum(<log(L)>d(beta)), while the second version of this is computed
    with half the number of integration points, and can be used to
    estimate the error in the evidence computation."""
    inds=np.argsort(betas)[::-1]

    sorted_logls=logls[inds]
    sorted_betas=betas[inds]

    dbetas=np.diff(np.concatenate((sorted_betas, [0.0])))
    dbetas2=np.diff(np.concatenate((sorted_betas[::2], [0.0])))

    return sorted_betas, -sorted_logls*dbetas, -sorted_logls[::2]*dbetas2

if __name__=='__main__':
    parser=ArgumentParser(description='Thermodynamic integration on PTMCMC samples.')
    
    parser.add_argument('files', metavar='OUTPUT_FILE', nargs='+', help='PTMCMC output files')
    parser.add_argument('--plotfile', metavar='FILE', default='evidence-integrand.pdf', help='plot output file')
    parser.add_argument('--evfile', metavar='FILE', default='evidence.dat', help='evidence output file')

    args=parser.parse_args()

    betas = np.array([1.0/extract_temp(f) for f in args.files])
    logls = np.array([get_mean_logl(f) for f in args.files])

    betas, integrands, integrands2 = thermo_integrands(logls, betas)

    evidence = np.sum(integrands)
    devidence = np.abs(evidence - np.sum(integrands2))

    pp.plot(betas, integrands/evidence)
    pp.xscale('log')
    pp.xlabel(r'$\beta$')
    pp.ylabel(r'$\frac{1}{\ln Z} \left \langle \log \mathcal{L} \right\rangle_\beta d\beta$')
    pp.savefig(args.plotfile)

    with open(args.evfile, 'w') as out:
        out.write('# ln-Z delta-ln-Z\n')
        out.write(str(evidence) + ' ' + str(devidence))
    
    
    

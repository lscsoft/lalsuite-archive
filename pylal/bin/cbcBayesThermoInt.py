#!/usr/bin/env python

from optparse import OptionParser
import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pp
import numpy as np
import scipy.integrate as si

matplotlib.rcParams['text.usetex'] = True

def skip_header(inp):
    skip = 0
    for line in inp:
        line = line.lstrip().split()
        skip += 1
        try:
            line.index('cycle')
            break
        except ValueError:
            pass

    return skip, line

def extract_temp(filename, burnin=0.5):
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
                temperature = inp.next().split()[i-1].strip()
                if temperature == 'variable':
                    variable = True
                    break
                else:
                    return float(temperature)
            except ValueError:
                pass

    with open(filename, 'r') as inp:
        if variable:
            skip, line = skip_header(inp)
            col = line.index('temp')
            data = np.loadtxt(filename, skiprows=skip, usecols=(col,))
            N = data.shape[0]

            return np.mean(data[N * burnin:])
        else:
            raise ValueError('extract_temp: did not find header line with \'Tchain\'')

def get_mean_logl(filename, burnin=0.5):
    """Returns the mean value of log(L) from the given filename,
    excluding the first 50% of samples as burn-in."""
    with open(filename, 'r') as inp:
        skip, line = skip_header(inp)

    col = line.index('logl')
    data = np.loadtxt(filename, skiprows=skip, usecols=(col,))
    N = data.shape[0]

    return np.mean(data[N * burnin:])

if __name__=='__main__':
    # Custom usage and help message
    usage = """%s [-h] [--plotfile FILE] [--evfile FILE] OUTPUT_FILE [OUTPUT_FILE ...]

Thermodynamic integration on PTMCMC samples.

positional arguments:
  OUTPUT_FILE      PTMCMC output files""" % (os.path.basename(sys.argv[0]))

    parser = OptionParser(usage=usage)
    parser.add_option('--plotfile', metavar='FILE', default='evidence-integrand.pdf', help='plot output file')
    parser.add_option('--evfile', metavar='FILE', default='evidence.dat', help='evidence output file')

    (options, args) = parser.parse_args()

    # Make positional arguments required
    if len(args)==0:
        parser.error('Positional filename argument(s) are required')

    betas = np.array([1.0/extract_temp(f) for f in args])
    logls = np.array([get_mean_logl(f) for f in args])

    inds = np.argsort(betas)[::-1]

    betas = betas[inds]
    logls = logls[inds]

    if betas[-1] != 0:
        betas2 = np.concatenate((betas[::2], [0]))
        betas = np.concatenate((betas, [0]))

        # Duplicate mean log-likelihood of hottest chain as a best guess for beta = 0.  This works
        # as long as the chains ahve extended to high enough temperature to sample the prior.
        logls2 = np.concatenate((logls[::2], [logls[-1]]))
        logls = np.concatenate((logls, [logls[-1]]))
    else:
        betas2 = np.concatenate((betas[:-1:2], betas[-1]))
        logls2 = np.concatenate((logls[:-1:2], logls[-1]))

    lnZ = -np.trapz(betas, logls)
    lnZ2 = -np.trapz(betas2, logls2)
    dlnZ = np.abs(lnZ - lnZ2)

    pp.plot(betas, betas*logls)
    pp.xscale('log')
    pp.xlabel(r'$\beta$')
    pp.ylabel(r'$\beta \left\langle \ln \mathcal{L} \right\rangle$')
    pp.savefig(options.plotfile)

    with open(options.evfile, 'w') as out:
        out.write('# ln-Z delta-ln-Z\n')
        out.write(str(lnZ) + ' ' + str(dlnZ))

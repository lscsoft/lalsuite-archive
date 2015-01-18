#!/usr/bin/env python

from optparse import OptionParser
import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pp
import numpy as np

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

def get_temperatures(filename):
    with open(filename, 'r') as file:
        skip, line = skip_header(file)
        cols = line.index('cycle'), line.index('temp')
        return np.loadtxt(filename, skiprows=skip, usecols=cols)

if __name__ == '__main__':
    usage = """%s [-h] [--plotfile FILE] OUTPUT_FILE [OUTPUT_FILE ...]

Plot evolution of dynamically adapted PTMCMC ladder.

positional arguments:
  OUTPUT_FILE      PTMCMC output files""" % (os.path.basename(sys.argv[0]))

    parser = OptionParser(usage=usage)
    parser.add_option('--plotfile', metavar='FILE', default='chain-evolution.pdf', help='Plot output file.')

    (options, args) = parser.parse_args()

    # Make positional arguments required
    if len(args)==0:
        parser.error('Positional filename argument(s) are required.')

    chains = sorted([get_temperatures(f) for f in args], key=lambda c: c[0, 0])
    for i, c in enumerate(chains[:4]):
        pp.plot(c[:,0], c[:, 1], label=r'$T_{%i}$' % i)
    pp.xlabel(r'No. of cycles')
    pp.ylabel(r'Temperature')
    pp.legend()
    pp.savefig(options.plotfile)

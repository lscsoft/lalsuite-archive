#!/usr/bin/env python

import itertools
import argparse
import re
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
    match = re.search(r'\.(\d+)$', filename)
    if match:
        index = int(match.group(1))
    else:
        index = None

    with open(filename, 'r') as file:
        skip, line = skip_header(file)
        cols = line.index('cycle'), line.index('temp')

        return index, np.loadtxt(filename, skiprows=skip, usecols=cols)

def get_ratios(filename, timeScale):
    match = re.search(r'\.(\d+)$', filename)
    if match:
        index = int(match.group(1))
    else:
        index = None

    with open(filename, 'r') as file:
        line = map(str.strip, file.next().split())
        cols = line.index('cycle'), line.index('swap_accepted')

        swaps = np.loadtxt(filename, skiprows=1, usecols=cols)
        bins = np.arange(np.min(swaps[:, 0]), np.max(swaps[:, 0]) + timeScale, timeScale)
        pairs = zip(np.digitize(swaps[:, 0], bins), swaps[:, 1])
        times = (bins[1:] + bins[:-1]) / 2.
        ratios = np.zeros(len(times))
        for i, g in itertools.groupby(pairs, lambda x: x[0]):
            g = list(x[1] for x in g)
            ratios[i - 1] = sum(g) / len(g)
        return index, times, ratios

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate diagnostic plot of chain temperature evolution.')
    parser.add_argument('outputFiles', nargs='+',
                        help='The chain output files.')
    parser.add_argument('--swap-files', nargs='*',
                        type=argparse.FileType('r'),
                        help='The swap statistic files for swap ratio plots.')
    parser.add_argument('--plot-file', type=argparse.FileType('w'),
                        help='The output plot file.',
                        default='chain-evolution.png')
    parser.add_argument('--indices', type=int,
                        nargs='*',
                        help='The indicies of the chains to plot.')
    parser.add_argument('--time-scale', type=int,
                        help='The averaging time-scale (in MCMC cycles) for acceptance ratios (default: 5000).',
                        default=5000)
    parser.add_argument('--log-time', action='store_true',
                        help='Plot time on a log-scale.')
    args = parser.parse_args()

    subplotCount = 2 if len(args.swap_files) > 0 else 1

    figure = pp.figure()
    axes = figure.add_subplot(subplotCount, 1, 1)

    chains = sorted(map(get_temperatures, args.outputFiles), key=lambda c: c[0])
    for i, c in chains:
        if args.indices is not None and i not in args.indices:
            continue

        axes.plot(c[:,0], c[:, 1], label=r'$T_{{{:}}}$'.format(i))
    axes.set_xlabel(r'No. of cycles')
    axes.set_ylabel(r'Temperature')
    axes.set_yscale('log')
    if args.indices is not None:
        axes.legend()
    xMin, xMax = axes.get_xlim()

    if len(args.swap_files) > 0:
        axes = figure.add_subplot(subplotCount, 1, 2)
        for file in args.swap_files:
            i, times, ratios = get_ratios(file.name, timeScale=args.time_scale)
            if args.indices is not None and i not in args.indices:
                continue

            axes.plot(times, ratios, label=r'$T_{{{:}}} \leftrightarrow T_{{{:}}}$'.format(i - 1, i))
        axes.set_xlabel(r'No. of cycles')
        axes.set_ylabel(r'Acceptance ratio')
        axes.set_xlim(xMin, xMax)
        axes.set_ylim(0, 1)
        if args.indices is not None:
            axes.legend()

    figure.savefig(args.plot_file.name)

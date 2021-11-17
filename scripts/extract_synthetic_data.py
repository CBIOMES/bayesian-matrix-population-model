#!/usr/bin/env python3

import os
import numpy as np
import netCDF4 as nc4
from shutil import copyfile
import argparse

default_templatedatafile = os.path.join(os.path.dirname(__file__), '..', 'data', 'size_distribution', 'template_synthetic_data.nc')

parser = argparse.ArgumentParser(description='Create a synthetic dataset from a sample simulation in a netCDF4 result file.')
parser.add_argument('resultfile', type=str, help='The name of the results file.')
parser.add_argument('index', type=int, help='The sample index used for synthetic data.')
parser.add_argument('syntheticdatafile', type=str, help='The name of output file.')
parser.add_argument('--model', type=str, default='m_bmb', help='The name of the model in the result file (default: "m_bmb").')
parser.add_argument('--templatedatafile', type=str, default=default_templatedatafile, help='A template data file with the correct values for PAR etc.; can be the original data file (default: "{}").'.format(default_templatedatafile))

args = parser.parse_args()

if os.path.isfile(args.syntheticdatafile):
    parser.error('"{}" already exists.'.format(args.syntheticdatafile))

copyfile(args.templatedatafile, args.syntheticdatafile)

with nc4.Dataset(args.syntheticdatafile, 'a') as ncd, nc4.Dataset(args.resultfile) as ncs:
    g = ncs[args.model]

    # new obs
    tmp = g['mod_obspos'][args.index, :, :].data
    ncd['w_obs'][:] = tmp / np.sum(tmp, axis=0)[None,:]

    timeindex = np.zeros(ncs['time'].size, dtype=bool)
    for t in ncs['obstime'][:].data:
        timeindex |= ncs['time'][:].data == t
    ncd['count'][:] = ncd['count'][0] * g['cell_count'][args.index, timeindex]
    ncd['abundance'][:] = ncd['abundance'][0] * g['cell_count'][args.index, timeindex]


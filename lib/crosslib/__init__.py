"""
Library for redshift-space density-momentum cross power spectrum
"""

import numpy as np
import h5py
import glob
from crosslib.util import _data_dir

import crosslib.power_spectrum

from crosslib.util import load_param
from crosslib.power_spectrum import load_linear_power, load_power_multipoles, compute_sigma_v

def load_lambda(isnp):
    filename = '%s/ps2d/ps2d_summary_%s.h5' % (_data_dir, isnp)

    d = {}
    summary = {}

    with h5py.File(filename, 'r') as f:
        d['k'] = f['k'][:]
        d['mu'] = f['mu'][:]
        d['lambda'] = f['lambda'][:]
        d['nrealisations'] = f['nrealisations'][()]

        for p in ['Pdd', 'Ppd']:
            summary[p] = f[p][:]
            summary['d' + p] = f['d' + p][:]

    d['summary'] = summary

    return d

"""
Library for redshift-space density-momentum cross power spectrum
"""

import numpy as np
import h5py
import glob
from crosslib.util import _data_dir

import crosslib.power_spectrum

def load_lambda(isnp):
    dirs = sorted(glob.glob('%s/ps2d/%s/0*' % (_data_dir, isnp)))
    if not dirs:
        raise FileNotFoundError('No subdirectories found in %s/%s/' %
                                (_data_dir, isnp))

    nrealisation = len(dirs)

    lambdas = np.arange(1, 101)*0.02
    nlambda = len(lambdas)
    
    Pdd = None
    Pdp = None
    
    for n, di in enumerate(dirs):
        for ilambda, lmbda in enumerate(lambdas):
            filename = '%s/ps2d_%.2f.h5' % (di, lmbda)
            with h5py.File(filename, 'r') as f:
                if Pdd is None:
                    shape = f['ps2d_dd'].shape
                    Pdd = np.empty((shape[0], shape[1], nlambda, nrealisation))
                    Pdp = np.empty_like(Pdd)
                    k = f['k'][:]
                    mu = f['mu'][:]

                Pdd[:, :, ilambda, n] = f['ps2d_dd'][:]
                Pdp[:, :, ilambda, n] = -f['ps2d_dp'][:]
                assert(abs(f['lambda'][()] - lmbda) < 1.0e-14)
    

    print(Pdd.shape)

    d = {}
    d['lambda'] = lambdas
    d['k'] = k
    d['mu'] = mu
    d['Pdd'] = Pdd
    d['Pdp'] = Pdp

    return d


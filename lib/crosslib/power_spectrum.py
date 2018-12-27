import numpy as np
import h5py
import glob

from crosslib.util import _data_dir

def load_multipoles(isnp):
    dirs = sorted(glob.glob('%s/ps/%s/0*' % (_data_dir, isnp)))
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
            filename = '%s/ps_%.2f.txt' % (di, lmbda)
            a = np.loadtxt(filename)

            if Pdd is None:
                Pdd = np.empty((a.shape[0], 3, nlambda, nrealisation))
                Pdp = np.empty_like(Pdd)
                k = a[:, 0]
                nmodes = np.zeros((a.shape[0]))

            Pdd[:, :, ilambda, n] = a[:, 1:4]
            Pdp[:, :, ilambda, n] = -a[:, 4:7]

        nmodes += a[:, 7]

    print(Pdd.shape)

    d = {}
    d['lambda'] = lambdas
    d['k'] = k
    d['nmodes'] = nmodes
    d['Pdd'] = Pdd
    d['Pdp'] = Pdp

    return d

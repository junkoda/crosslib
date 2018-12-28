import numpy as np
import math
import h5py
import glob

from crosslib.util import _data_dir

def load_multipoles(isnp):
    """
    Returns auto power Pdd and cross power Ppd as a function of lambda

    Args:
      isnp (str): snapshot index 000 - 010

    Returns:
      d (dict):
        d['lambda'] (array): lambda[ilambda]
        d['k'] (array): k[ik, imu] [h/Mpc]
        d['nmodes'] (array): number of independent k modes (sum of n realisations)
        d['Pdd'] (array): Pdd[ik, imu, ilambda, irealisation]
        d['Ppd'] (array): Ppd[ik, imu, ilambda, irealisation]
        d['nrealisations']: number of realisations n
        
        d['summary'] (dict): mean and error in the mean of n realisations
        d['summary']['Pdd'] (array):  mean Pdd[ik, imu, ilambda]
        d['summary']['dPdd'] (array): error in the mean dPdd[ik, imu, ilambda]
        d['summary']['Ppd'] (array):  mean Pdd[ik, imu, ilambda]
        d['summary']['dPpd'] (array): error in the mean dPpd[ik, imu, ilambda]

    Notes:
      k is mean wavenumber in the bin
    """
    dirs = sorted(glob.glob('%s/ps/%s/0*' % (_data_dir, isnp)))
    if not dirs:
        raise FileNotFoundError('No subdirectories found in %s/%s/' %
                                (_data_dir, isnp))

    nrealisations = len(dirs)

    lambdas = np.arange(1, 101)*0.02
    nlambda = len(lambdas)
    
    Pdd = None
    Ppd = None
    
    for n, di in enumerate(dirs):
        for ilambda, lmbda in enumerate(lambdas):
            filename = '%s/ps_%.2f.txt' % (di, lmbda)
            a = np.loadtxt(filename)

            if Pdd is None:
                Pdd = np.empty((a.shape[0], 3, nlambda, nrealisations))
                Ppd = np.empty_like(Pdd)
                k = a[:, 0]
                nmodes = np.zeros((a.shape[0]))

            Pdd[:, :, ilambda, n] = a[:, 1:4]
            Ppd[:, :, ilambda, n] = -a[:, 4:7]

        nmodes += a[:, 7]

    print(Pdd.shape)

    d = {}
    d['lambda'] = lambdas
    d['k'] = k
    d['nmodes'] = nmodes
    d['Pdd'] = Pdd
    d['Ppd'] = Ppd
    d['nrealisations'] = nrealisations

    summary = {}
    for p in ['Pdd', 'Ppd']:
        summary[p] = np.mean(d[p], axis=3)
        summary['d' + p] = np.std(d[p], axis=3)/math.sqrt(nrealisations)

    d['summary'] = summary
    
    return d

def load_lambda_all(isnp):
    dirs = sorted(glob.glob('%s/ps2d/%s/0*' % (_data_dir, isnp)))
    if not dirs:
        raise FileNotFoundError('No subdirectories found in %s/%s/' %
                                (_data_dir, isnp))

    nrealisation = len(dirs)

    lambdas = np.arange(1, 101)*0.02
    nlambda = len(lambdas)
    
    Pdd = None
    Ppd = None
    
    for n, di in enumerate(dirs):
        for ilambda, lmbda in enumerate(lambdas):
            filename = '%s/ps2d_%.2f.h5' % (di, lmbda)
            with h5py.File(filename, 'r') as f:
                if Pdd is None:
                    shape = f['ps2d_dd'].shape
                    Pdd = np.empty((shape[0], shape[1], nlambda, nrealisation))
                    Ppd = np.empty_like(Pdd)
                    k = f['k'][:]
                    mu = f['mu'][:]

                Pdd[:, :, ilambda, n] = f['ps2d_dd'][:]
                Ppd[:, :, ilambda, n] = -f['ps2d_dp'][:]
                assert(abs(f['lambda'][()] - lmbda) < 1.0e-14)
    

    print(Pdd.shape)

    d = {}
    d['lambda'] = lambdas
    d['k'] = k
    d['mu'] = mu
    d['Pdd'] = Pdd
    d['Ppd'] = Ppd

    return d


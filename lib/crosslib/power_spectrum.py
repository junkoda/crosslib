import numpy as np
import math
import h5py
import glob
import numbers
     
from crosslib.util import _data_dir, load_param


def load_linear_power(isnp=None, *, k_min=None, k_max=None):
    """
    Linearly extrapolated power spectrum at redshift of isnp
    Returns z=0 power spectrum if isnp = None

    Args:
      isnp (str):    snapshot index
      k_min (float): minimum k
      k_max (float): maximum k

    Returns:
      d (dict)
        d['k'] (array): wavenumbeer [h/Mpc]
        d['P'] (array): linear power spectrum [1/h Mpc]^3
        d['interp'] (function): interpolation function P(k)
    """

    filename = '%s/planck_matterpower.dat' % (_data_dir,)
    a = np.loadtxt(filename)

    if isnp is not None:
        param = load_param(isnp)
        a[:, 1] *= param['D']**2

    d = {}
    d['a'] = a

    # reduce k range to [k_min, k_max]
    if k_min is not None or k_max is not None:
        if k_min is None:
            k_min = 0.0
        if k_max is None:
            k_max = math.inf

        idx = np.logical_and(k_min <= a[:, 0], a[:, 0] <= k_max)
        d['k'] = a[idx, 0]
        d['P'] = a[idx, 1]
    else:
        d['k'] = a[:, 0]
        d['P'] = a[:, 1]

    # interpolated function
    try:
        from scipy.interpolate import interp1d
        d['interp'] = interp1d(a[:, 0], a[:, 1], kind='cubic')
    except ImportError:
        print('Warning: scipy.interpolate unabailable '
              'for linear power interpolation.')

    return d


def load_power_multipoles(isnp):
    """
    Returns auto power Pdd and cross power Ppd as a function of lambda

    Args:
      isnp (str): snapshot index 000 - 010

    Returns:
      d (dict):
        d['lambda'] (array): lambda[ilambda]
        d['k'] (array): k[ik, imu] [h/Mpc]
        d['nmodes'] (array): number of independent k modes
                             (sum of n realisations)
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

    lambdas = np.arange(0, 101)*0.02
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

def load_power2d(isnp, *, lmbda=1.00, unit='displacement'):
    """
    Load 2D power spectra for one lambda

    Args:
      lmbda (float): lambda value
      unit (str): unit of the velocity 'km/s' or 'displacement' (1/hMpc)

    Returns:
      d (dict)
       d['k']      (array): k[ik, imu]
       d['mu']     (array): mu[ik, imu]
       d['lambda'] (float)
       d['Pdd'], d['Pdp'], d['Ppp'] (array): Pab[ik, imu, irealisation]
       d['summary'] (dict):
         Pdd, Ppd, Ppp    (array): Pab[ik, imu] mean
         dPdd, dPpd, dPpp (array): standard deviation in the mean
    Note:
      Use load_lambda() for all lambdas
    """

    if not isinstance(lmbda, str):
        lmbda = '%.2f' % lmbda
    
    data_dir = '/Users/junkoda/Research/cross/doraemon/cross1'
    
    filenames = glob.glob('%s/ps2d/%s/0*/ps2d_%s.h5' % (data_dir, isnp, lmbda))
    filenames = sorted(filenames)

    if not filenames:
        raise FileNotFoundError('No files found in %s/ps2d/%s/' % (data_dir, isnp))

    n = len(filenames)
    print(n, 'realisations')

    # Returning data
    d = {}
    summary = {}

    d['nrealisations'] = n
    d['summary'] = summary
    
    P2_dd = None
    P2_pd = None
    #P2_pp = None

    for i, filename in enumerate(filenames):
        with h5py.File(filename, 'r') as f:
            if P2_dd is None:
                shape = f['ps2d_dd'].shape
                P2_dd = np.empty((shape[0], shape[1], n))
                P2_pd = np.empty_like(P2_dd)
                #P2_pp = np.empty_like(P2_dd)
                
                d['k'] = f['k'][:]
                d['mu'] = f['mu'][:]
                d['lambda'] = f['lambda'][()]

            P2_dd[:, :, i] = f['ps2d_dd'][:]
            P2_pd[:, :, i] = -f['ps2d_dp'][:]
            #P2_pp[:, :, i] = f['ps2d_pp'][:]

    d['Pdd'] = P2_dd
    d['Ppd'] = P2_pd

    if unit == 'displacement':
        pass
    elif unit == 'km/s':
        param = load_param(isnp)
        aH = param['a']*param['H']
        d['Ppd'] *= aH
    else:
        raise ValueError('Unknown unit: %s' % unit)
    
    #d['Ppp'] = fac_vel**2*P2_pp

    for p in ['Pdd', 'Ppd']:
        summary[p] = np.mean(d[p], axis=2)
        if n > 0:
            summary['d' + p] = np.std(d[p], axis=2)/math.sqrt(n)

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

    d = {}
    d['lambda'] = lambdas
    d['k'] = k
    d['mu'] = mu
    d['Pdd'] = Pdd
    d['Ppd'] = Ppd

    return d


def compute_sigma_v(isnp):
    """
    compute linar sigma_v = \int P(k) dk/(6 pi^2)
    """
    linear = load_linear_power(None)
    
    param = load_param(isnp)
    fac = param['f']*param['D']

    k = linear['k']
    P = linear['P']
    dk = k[1:] - k[:-1]

    # Trapezoidal integral
    sigma2v = 0.5*np.sum((P[1:] + P[:-1])*(k[1:] - k[:-1]))/(6.0*math.pi**2)

    return fac*math.sqrt(sigma2v)

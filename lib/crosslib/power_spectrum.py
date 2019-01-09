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

def load_power2d(isnp, lmbda):
    """
    Load 2D power spectra for one lambda

    Args:
      lmbda (float): lambda value
                     0 for real space, 1 for redshift space

    Returns:
      d (dict)
       d['k']      (array): k[ik, imu]
       d['mu']     (array): mu[ik, imu]
       d['lambda'] (float)
       d['Pdd'], d['Pdp'], d['Ppp'] (array): Pab[ik, imu, irealisation]
       d['summary'] (dict):
         Pdd, Ppd, Ppp    (array): Pab[ik, imu] mean
         dPdd, dPpd, dPpp (array): standard deviation in the mean
    """

    if not isinstance(lmbda, str):
        lmbda = '%.2f' % lmbda
    
    filenames = glob.glob('%s/ps2d/%s/0*/ps2d_%s.h5' % (_data_dir, isnp, lmbda))
    filenames = sorted(filenames)

    if not filenames:
        raise FileNotFoundError('No ps2d_*.h5 files found in %s/ps2d/%s/'
                                % (_data_dir, isnp))

    n = len(filenames)

    # Returning data
    d = {}
    summary = {}

    d['nrealisations'] = n
    d['summary'] = summary
    
    P2_dd = None
    P2_pd = None
    P2_pp = None

    for i, filename in enumerate(filenames):
        with h5py.File(filename, 'r') as f:
            if P2_dd is None:
                shape = f['ps2d_dd'].shape
                P2_dd = np.empty((shape[0], shape[1], n))
                P2_pd = np.empty_like(P2_dd)
                P2_pp = np.empty_like(P2_dd)
                
                d['k'] = f['k'][:]
                d['mu'] = f['mu'][:]
                d['lambda'] = f['lambda'][()]

            P2_dd[:, :, i] = f['ps2d_dd'][:]
            P2_pd[:, :, i] = -f['ps2d_dp'][:] # Ppd = - Pdp
            P2_pp[:, :, i] = f['ps2d_pp'][:]

    d['Pdd'] = P2_dd
    d['Ppd'] = P2_pd
    d['Ppp'] = P2_pp

    for p in ['Pdd', 'Ppd', 'Ppp']:
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


def load_halofit_power(isnp):
    """
    Args:
      isnp (str): snasphot index
    """
    filename = '%s/power_spectrum/%s/halofit_matterpower.dat' % (_data_dir, isnp)
    a = np.loadtxt(filename)
    
    d = {}
    d['k'] = a[:, 0]
    d['P'] = a[:, 1]

    try:
        from scipy.interpolate import interp1d
        d['interp'] = interp1d(d['k'], d['P'], kind='cubic')
    except ImportError:
        print('Warning: scipy.interpolate unabailable '
              'for linear power interpolation.')


    return d


def load_theta_power_bel(isnp, *, Pdd=None, linear=None,
                         Ptt_simple=False):
                                
    """
    Density- Velocity-divergence cross power using Bell et al. formula
    https://arxiv.org/abs/1809.09338

    Args:
      sim (str): simulation name
      isnp (str): snapshot index
      Pdd (dict): [optional] Non-linear Pdd(k, z) dictionary with 'k' and 'P'
                  halofit is loaded if not provided
      linear (dict): [optional] linear P(k, z) dictionary with 'k' and 'P'
                  linear is loaded if not provided

    Returns:
      d['k']: wavenumber kdd if provided, linear['k'] otherwise
      d['Pdd']: Pdd if provided, linear['P'] otherwise
      d['Pdt']: P_delta_theta
      d['Ptt']: P_theta_theta

    Reference:
      Bell et al. https://arxiv.org/abs/1809.09338
    """
    
    param = load_param()
    D = param['snapshot'][isnp]['D']
    sigma8 = D*param['sigma_8']

    a1 = -0.817 + 3.198*sigma8
    a2 = 0.877 - 4.191*sigma8
    a3 = -1.199 + 4.629*sigma8
    kd_inv = -0.111 + 3.811*sigma8**2
    b = 0.091 + 0.702*sigma8
    kt_inv = -0.048 + 1.917*sigma8**2

    if linear is None:
        linear = load_linear_power(isnp)

    if Pdd is None:
        Pdd = load_halofit_power(isnp)
    
    if not ('k' in linear and 'P' in linear):
        raise ValueError("dict linear must contain 'k' and 'P'")

    if not ('k' in Pdd and 'P' in Pdd):
        raise ValueError("dict Pdd must contain 'k' and 'P'")

    k = Pdd['k']
    Pdd = Pdd['P']

    # Have P and Pdd at same k
    if not np.all(k == linear['k']):
        from scipy.interpolate import interp1d
        P = interp1d(linear['k'], linear['P'], kind='cubic')(k)
    else:
        P = linear['P']
    
    d = {}
    d['k'] = k
    d['Pdd'] = Pdd
    d['Pdt'] = np.sqrt(Pdd*P)*np.exp(-k*kd_inv - b*k**6)
    
    if Ptt_simple:
        d['Ptt'] = P*np.exp(-k*kt_inv)
    else:
        d['Ptt'] = P*np.exp(-k*(a1 + a2*k + a3*k**2))

    return d

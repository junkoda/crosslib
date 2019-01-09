import os
import json

try:
    # The data directory can be defined in crosslib_dir.py
    import crosslib.crosslib_dir
    _cross_dir = crosslib.crosslib_dir._cross_dir
    _data_dir = _cross_dir + '/data'
except ImportError:
    # if crosslib_dir.py is not provided, assume data are in ../../data/
    _this_dir = os.path.dirname(os.path.realpath(__file__))
    _data_dir = os.path.abspath(_this_dir + '/../../data')


def load_param(isnp=None):
    """
    Return simulation info
    """

    with open('%s/param.json' % (_data_dir)) as f:
        d = json.load(f)

    if isnp is None:
        return d

    if isnp not in d['snapshot']:
        raise ValueError('isnp %s is not available.' % (isnp))

    return d['snapshot'][isnp]

import os

try:
    # The data directory can be defined in crosslib_dir.py
    import crosslib.crosslib_dir
    _data_dir = crosslib.crosslib_dir._data_dir
except ImportError:
    # if crosslib_dir.py is not provided, assume data are in ../../data/
    _this_dir = os.path.dirname(os.path.realpath(__file__))
    _data_dir = os.path.abspath(_this_dir + '/../../data')
    
#_data_dir = '/Users/junkoda/Research/cross/doraemon/cross'



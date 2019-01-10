#
# usage
# $ pip3 install .
#

from setuptools import setup
import os

this_dir = os.path.dirname(os.path.realpath(__file__))
cross_dir = os.path.abspath(this_dir + '/../..')

with open('crosslib/cross_dir.py', 'w') as f:
    f.write("_cross_dir='%s'\n" % cross_dir)

setup(
    name='crosslib',
    version='0.0.1',
    python_requires='>=3',
    install_requires=['numpy', 'scipy', 'h5py'],
    author='Jun Koda',
    url='https://github.com/junkoda/crosslib'
)


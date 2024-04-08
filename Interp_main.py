# Copyright (c) Xie Zhang
# Distributed under the terms of the MIT License.

import numpy as np
from transmatrix_interp.Transmatrix import Transmatrix
import argparse

''' Interpolation module for eigenvalue and transition matrix

This module is the main function for interpolation.
different interpolation schemes like cubic, linear 
and wannier interpolation are included. The wannier
interpolation are executed through VASP
'''

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-p', '--path', type=str, default='./test_files/', help='Path of transmatrix file, should include IBZKPT file')
    parser.add_argument('-VB', '--nvb', type=int, nargs='+', default=[32, 33, 34, 35, 36], help='Number of valence bands')
    parser.add_argument('-CB', '--ncb', type=int, nargs='+', default=[37, 38, 39, 40, 41],help='Number of conduction bands')
    parser.add_argument('-m', '--magnification', type=int, nargs=3, default=[10, 10, 10],help='Mesh to generate kpoints (nx, ny, nz)')
    parser.add_argument('-M', '--method', nargs='+', choices=['linear', 'cubic'], default=['linear', 'cubic'], help='Interpolate methods first for trans second for eigenvalues')
    args = parser.parse_args()

    a = np.array(args.nvb)
    b = np.array(args.ncb)
    root = args.path
    mag = args.magnification
    method = args.method

    interp = Transmatrix(a, b, root, method)
    interp.write_to_file(mag[0], mag[1], mag[2])

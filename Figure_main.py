# Copyright (c) Xie Zhang
# Distributed under the terms of the MIT License.

import numpy as np
from rate.figure_plot import figure_plot
import argparse

''' Visualization module for all calculation results
This module is the main function for visualization of
previous calculation results. The visualization includes 
the absorption coefficients, radiative recombination 
coefficients and carrier distribution
'''
if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-n', '--name', type=str, default='GaAs', help='the name of materials')
    parser.add_argument('-V', '--Vcell', type=float, nargs=1, default=300.0, help='volume of primitive cell, found in outcar')
    parser.add_argument('-S', '--Soc', type=bool,default=False, help='whether soc')
    parser.add_argument('-p', '--path', type=str, default='./test_files/', help='dimension of materials')
    parser.add_argument('-T', '--Trange', type=float, nargs='+', default=[300], help='range of temperature')
    parser.add_argument('-D', '--Density', type=float, nargs='+', default=[1E19], help='range of density')
    parser.add_argument('-VB', '--Valence', type=int, nargs='+', default=[32, 33, 34, 35, 36], help='valence bands')
    parser.add_argument('-CB', '--Conduction', type=int, nargs='+', default=[37, 38, 39, 40, 41], help='conduction bands')
    parser.add_argument('-d', '--direction', type=str, default='z', help='direction of chosen plane')
    parser.add_argument('-k', '--kpoints', type=float, nargs='+', default=[0.0], help='kpoint coordiante')
    parser.add_argument('-E', '--Eg', type=float, nargs=1, default=2.39, help='Auxiliary line')
    parser.add_argument('-B', '--startband', type=int, nargs=1, default = 20, help='the start num of band involving calculation, '
                                                                                  'used for wannier interpolation')

    args = parser.parse_args()

    Name = args.name
    path = args.path
    Vcell = args.Vcell
    T = args.Trange
    D = args.Density
    VB = args.Valence
    CB = args.Conduction
    soc = args.Soc
    direction = args.direction
    kpoints = args.kpoints
    Eg = args.Eg
    start_band = args.startband

    Figure = figure_plot(Name, path)
    Figure.adsorption_coeff_plot(Eg= Eg)
    #Figure.radiative_coeff_plot(1E15, 1E19)
    #Figure.plot_k_slice(direction, kpoints, Vcell, T, D, VB, CB, soc, start_band)
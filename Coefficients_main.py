# Copyright (c) Xie Zhang
# Distributed under the terms of the MIT License.

import numpy as np
from rate.rate_calculation import rate_cal
import argparse

''' Core module for coefficients calculation
This module is the main function for absorption coefficients
and radiative recombination coefficients(B coefficiens) calculation
'''

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-r', '--refractive', type=float, default=1.58, help=' refractive index !!!!')
    parser.add_argument('-n', '--name', type=str, default='GaAs', help='the name of materials')
    parser.add_argument('-V', '--Vcell', type=float, nargs=1, default=183.0, help='volume of primitive cell, found in outcar')
    parser.add_argument('-S', '--Soc', type=bool,default=True, help='whether soc')
    parser.add_argument('-p', '--path', type=str, default='./test_files/', help='the work path')
    parser.add_argument('-T', '--Trange', type=float, nargs='+', default=[100, 200, 300], help='range of temperature')
    parser.add_argument('-D', '--Density', type=float, nargs='+', default=[1E15, 1E16, 1E17, 1E18, 1E19], help='range of density')
    parser.add_argument('-sig', '--sigma', type=float, nargs=1, default=0.001, help='the sigma coefficients for gauss function')
    parser.add_argument('-E', '--Energy', type=float, nargs='+', default=[1, 3, 301],
                        help='range of Energy for absorption calculation')
    parser.add_argument('-B', '--startband', type=int, nargs=1, default = 20, help='the start num of band involving calculation, '
                                                                                  'used for wannier interpolation')

    args = parser.parse_args()

    mu = args.refractive
    Vcell = args.Vcell
    soc = args.Soc
    name = args.name
    Erange = args.Energy
    start_band = args.startband

    root = args.path
    Density = np.array(args.Density)
    Trange = np.array(args.Trange)

    Rate = rate_cal(name, mu, Vcell, False, root, start_band)
    Rate.absorption(Erange[0], Erange[1], Erange[2], 0.024)
    #Rate.get_Radcoeffi_results(Trange, Density)
    #Rate.radiative_coeff_plot(Trange, Density)

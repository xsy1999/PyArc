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

    parser.add_argument('-A', '--absorp', type=lambda x: (str(x).lower() == 'true'), default=True, help='whether soc')
    parser.add_argument('-R', '--radiative', type=lambda x: (str(x).lower() == 'true'), default=True, help='whether soc')
    parser.add_argument('-u', '--refractive', type=float, default=3.9476, help=' refractive index !!!!')
    parser.add_argument('-n', '--name', type=str, default='GaAs', help='the name of materials')
    parser.add_argument('-V', '--Vcell', type=float, default=183.0, help='volume of primitive cell, found in outcar')
    parser.add_argument('-S', '--Soc', type=lambda x: (str(x).lower() == 'true'), default=False, help='whether soc')
    parser.add_argument('-p', '--path', type=str, default='./test_files/', help='the work path')
    parser.add_argument('-T', '--Trange', type=float, nargs='+', default=[100, 200, 300], help='range of temperature')
    parser.add_argument('-D', '--Density', type=float, nargs='+', default=[1E15, 1E16, 1E17, 1E18, 1E19], help='range of density')
    parser.add_argument('-sig', '--sigma', type=float, default=0.024, help='the sigma coefficients for gauss function')
    parser.add_argument('-E', '--Energy', type=float, nargs='+', default=[1, 3, 301],
                        help='range of Energy for absorption calculation, the last number is the total sampling steps in energy range')
    parser.add_argument('-B', '--startband', type=int, default=20, help='the start num of band involving calculation, '
                                                                                  'used for wannier interpolation')

    args = parser.parse_args()

    mu = args.refractive
    Vcell = args.Vcell
    soc = args.Soc
    name = args.name
    Erange = args.Energy
    start_band = args.startband
    absorp = args.absorp
    radiative = args.radiative
    sigma = args.sigma

    root = args.path
    Density = np.array(args.Density)
    Trange = np.array(args.Trange)

    Rate = rate_cal(name, mu, Vcell, soc, root, start_band)

    if absorp == True:
        Rate.absorption(Erange[0], Erange[1], Erange[2], sigma)
    if radiative == True:
        Rate.get_Radcoeffi_results(Trange, Density)
    #Rate.radiative_coeff_plot(Trange, Density)

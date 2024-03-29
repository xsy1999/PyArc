# Copyright (c) Xie Zhang
# Distributed under the terms of the MIT License.

import numpy as np
import os
import re

''' Input module

Input module offers several functions that help to 
read different files generated by VASP.
'''

#code for txt file input

def read_EIGEN(file, soc=False) -> (float, np.ndarray, np.ndarray):
    if os.path.exists(file) is not True:
        raise IOError("The %s file does not exist" %file)
    filedata = open(file)
    EIG = filedata.readlines()
    Info = EIG[5].split()
    Info = list(map(int, Info))
    if soc is True:
        vb = Info[0]
    else:
        vb = int(Info[0]/2)
    bands = Info[2]
    kpts = Info[1]
    data = []
    for num, line in enumerate(EIG):
        eachline = line.split()
        if num > 6 and line != ' \n':
            value = list(map(float, eachline[0 : 3]))
            data.append(value)
    data = np.array(data)
    Eigenvalue = data.reshape(kpts, bands+1, 3)
    order = np.lexsort((Eigenvalue[:, 0, 2], Eigenvalue[:, 0, 1], Eigenvalue[:, 0, 0]))
    Eigenvalue = Eigenvalue[order]
    kpoints = np.round(Eigenvalue[:, 0, :], 5)
    ksize = [np.unique(kpoints[:, 0]).shape[0], np.unique(kpoints[:, 1]).shape[0], np.unique(kpoints[:, 2]).shape[0]]
    Eigenvalue = Eigenvalue[:, 1::, 1]
    evbm = max(Eigenvalue[:, vb-1])
    ksize.append(bands)
    Eigenvalue = Eigenvalue.reshape(ksize)
    return ksize, evbm, vb, Eigenvalue

def read_kpt(file) -> (int, np.ndarray):
    if os.path.exists(file) is not True:
        raise IOError("The %s file does not exist" %file)
    Kptfile = open(file)
    Kpt = Kptfile.readlines()[3::]
    num = len(Kpt)
    kpoints = np.array([list(map(float, Kpt[i].split())) for i in range(num)])
    return num, kpoints

def read_transmatrix(file) -> (int, list):
    if os.path.exists(file) is not True:
        raise IOError("The %s file does not exist" %file)
    Transfile = open(file)
    Transdata = Transfile.readlines()
    num = len(Transdata)
    Trans = [list(map(float,Transdata[i].split())) for i in range(num)]
    return num, Trans

def read_radiative(file) -> ():
    if os.path.exists(file) is not True:
        raise IOError("The %s file does not exist" % file)
    density = []
    Bcoeff = []
    lifetimes = []
    Bcoef_file = open(file)
    Bcoef_data = Bcoef_file.readlines()
    for i in range(len(Bcoef_data)):
        if i == 0:
            pattern = r'\[(.*?)\]'
            matches = re.findall(pattern, Bcoef_data[i])[0].split()
            T = np.sort(np.array([float(char) for char in matches]))
        if i > 0:
            data = Bcoef_data[i].replace('[', '').replace(']', '').split()
            results = [list(map(float, data[j].split())) for j in range(len(data))]
            density.append(results[0])
            Bcoeff.append(results[1: int(T.shape[0]) + 1])
            lifetimes.append(results[1 + int(T.shape[0]): 2 * int(T.shape[0]) + 1])
    density = np.squeeze(np.array(density))
    Bcoeff = np.squeeze(np.array(Bcoeff)).transpose(1, 0)
    lifetimes = np.squeeze(np.array(lifetimes)).transpose(1, 0)
    return T, density, Bcoeff, lifetimes

def read_adsorption(file):
    if os.path.exists(file) is not True:
        raise IOError("The %s file does not exist" % file)
    adsorp_file = open(file)
    adsorp_data = adsorp_file.readlines()[1::]
    adsorp_coeff = [list(adsorp_data[i].replace('[', '').replace(']', '').split())\
                    for i in range(len(adsorp_data))]
    adsorp_coeff = np.array([list(map(float, sublist)) for sublist in adsorp_coeff])
    Erange = adsorp_coeff[:, 0]
    alpha = adsorp_coeff[:, 1::]
    return Erange, alpha


def read_wannier_eigs(file) -> (list, int, np.ndarray):
    if os.path.exists(file) is not True:
        raise IOError("The %s file does not exist" %file)
    eigtest = np.loadtxt(file)
    Eigfile = open(file)
    Eig = Eigfile.readlines()
    Eigdim = Eig[1].split()[-3::]
    Eigdim = [int(Eigdim[i]) for i in range(len(Eigdim))]
    num = eigtest.shape[0]
    Eigs = eigtest.reshape(Eigdim[0] * Eigdim[1] * Eigdim[2], -1, 5)
    return Eigdim, num, Eigs

# Copyright (c) Xie Zhang
# Distributed under the terms of the MIT License.

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
from typing import Optional, Union
import mpl_toolkits.axisartist.floating_axes as floating_axes
from transmatrix_interp.interpolator import bicubic_interp, linear_interp
from mpl_toolkits.axisartist.grid_finder import (FixedLocator, DictFormatter)
from rate.constants import const
from input import read_EIGEN, read_kpt, read_transmatrix, read_wannier_eigs
import os
from tqdm import tqdm

'''Calculation module for absorption and radiative coefficients

In this module, all the parameters and eigenvalues and transition
matrix values are combined to obtain the final absorption and 
radiative recombination coefficients.
'''

class rate_cal:
    def __init__(self, name: str, mu : float, Vcell : float, soc=False, root= './', start_band=0):
        self.path = root
        self.name = name
        self.mu = mu
        self.Vcell_cm = Vcell * const.Ang2cm ** 3
        self.Vcell_Bohr = Vcell * (const.AngtoBohr ** 3)
        self.conversion_factor = const.eV2erg * (const.hbar / const.abohr) ** 2
        FineKptfile = os.path.join(self.path + 'wannier90_geninterp.kpt')
        FineEigfile = os.path.join(self.path + 'wannier90_geninterp.dat')
        FineTransfile = os.path.join(self.path + 'matrix_fine.dat')
        Ori_eigfile = os.path.join(self.path + 'EIGENVAL')
        self.soc = soc
        self.nktot, kpts = read_kpt(FineKptfile)
        ksize, _, self.vborder, Eigenvalue = read_EIGEN(Ori_eigfile, self.soc)
        self.interpolator = bicubic_interp(1)
        self.kmesh = np.array([len(np.unique(kpts[:, 1])), len(np.unique(kpts[:, 2])),\
                      len(np.unique(kpts[:, 3]))])

        #substitute the kpoints coordinates in wannier.dat
        #obtain dense mesh Eigenvalues
        if os.path.exists(FineEigfile) is True:
            _, eignum, self.Eigs = read_wannier_eigs(FineEigfile)
        elif os.path.exists(Ori_eigfile):
            Eigs = self.interpolator.interpolate_kpoints(Eigenvalue, kpts[:, 1:])
            Eigs = np.expand_dims(Eigs, axis=2)
            eignum = Eigs.shape[0] * Eigs.shape[1]
            smear = np.zeros((Eigs.shape[0], Eigs.shape[1], 4))
            self.Eigs = np.concatenate((smear, Eigs), axis=2)
        else:
            raise IOError("no Eigenvalue file of dense mesh")

        for i in range(self.nktot):
            self.Eigs[i, :, 1:4] = kpts[i, 1:4]

        Eigs = self.Eigs[:, :, -1]
        self.evbm = max(Eigs[:, self.vborder-start_band-1])
        self.ecbm = min(Eigs[:, self.vborder-start_band])
        self.Eigs[:, :, -1] -= self.evbm

        self.nb = int(eignum / self.nktot)
        self.num, Trans = read_transmatrix(FineTransfile)

        #write all information into one file transmatrix.dat
        transtotal = []
        if Path(self.path + 'transmatrix_new.dat').is_file():
            file = self.path + 'transmatrix_new.dat'
            transtotal = np.loadtxt(file)
        else:
            f = open(self.path + 'transmatrix_new.dat', 'w')
            with tqdm(total=len(Trans)) as pbar:
                for num, line in enumerate(Trans):
                    pbar.set_description('Combining data:')
                    data = []
                    for i in range(len(line)):
                        if i < 2:
                            data.append(float(line[i]))
                        elif i ==2:
                            data.append(round(float(Eigs[num % self.nktot, int(line[0]) - start_band - 1]) - self.evbm, 10))
                            data.append(round(float(Eigs[num % self.nktot, int(line[1]) - start_band - 1]) - self.evbm, 10))
                            data.append(float(line[i]))
                        else :
                            data.append(float(line[i]))

                    transtotal.append(data)

                    f.write("%s \n" % str(data).replace(',', ' ').replace('[', '').replace(']', ' '))

                    pbar.update(1)
            f.close()
            transtotal = np.array(transtotal)
        self.total = transtotal
        self.matrix = np.array(transtotal[:, -3::])
        self.ebc = np.array(transtotal[:, 2]) # the energy of bands combination, ebc for the top one
        self.ebv = np.array(transtotal[:, 3]) # the energy of bands combination, ebv for the bottom one
        self.cb = list(set(transtotal[:, 0]))
        self.vb = list(set(transtotal[:, 1]))
        self.ncb = len(self.cb)
        self.nvb = len(self.vb)

    def deltagauss(self, dE, sigma):
        ans = np.exp(-dE ** 2 / sigma ** 2) / np.sqrt(np.pi) / sigma
        return ans

    #calculation of absorption coefficients
    def absorption(self, Emin, Emax, Enum, sigma: Optional[float] = None):
        #check whether the transmitrix.dat file exist
        if Path(self.path + 'transmatrix.dat').is_file():
            file = self.path + 'transmatrix.dat'
            transtotal = np.loadtxt(file)
            matrix = np.array(transtotal[:, -3::])
            kw = 1 / transtotal.shape[0] * len(set(matrix[:, 3])) * len(set(matrix[:, 4])) #!!!how to set
            ebv = np.array(transtotal[:, 3])
            ebc = np.array(transtotal[:, 2])
        else:
            ebv = self.ebv
            ebc = self.ebc
            matrix = self.matrix
            kw = 1 / self.nktot
        
        #self.Vcell_Bohr = 5.9698164876 ** 3 * (const.AngtoBohr ** 3)
        n_r = self.mu

        prefactor = 4 * np.pi ** 2 / self.Vcell_Bohr

        Erange = np.linspace(Emin, Emax, num=int(Enum))

        alpha =[]

        epsilon_ave = 0
        for en in tqdm(Erange, desc="Erange"):
            # n_r = np.sqrt(7.871 * 3.029 / (3.029**2 - en**2) + 1) # C. He et al., Cryst. Res. Technol. 2019, 54, 1900011
            freq = en * const.eVtoHar
            deltaE = (ebc - ebv - en) * const.eVtoHar
            epsilon_x_tmp = kw * matrix[:, 0] ** 2 * self.deltagauss(deltaE, sigma * const.eVtoHar)
            epsilon_y_tmp = kw * matrix[:, 1] ** 2 * self.deltagauss(deltaE, sigma * const.eVtoHar)
            epsilon_z_tmp = kw * matrix[:, 2] ** 2 * self.deltagauss(deltaE, sigma * const.eVtoHar)
            epsilon_ave_tmp = kw * (matrix[:, 0] ** 2 + matrix[:, 1] ** 2 + matrix[:, 2] ** 2) \
                    / 3.0 * self.deltagauss(deltaE, sigma * const.eVtoHar)
            epsilon_x = prefactor / freq ** 2 * epsilon_x_tmp.sum()
            epsilon_y = prefactor / freq ** 2 * epsilon_y_tmp.sum()
            epsilon_z = prefactor / freq ** 2 * epsilon_z_tmp.sum()
            epsilon_ave = prefactor / freq ** 2 * epsilon_ave_tmp.sum()
            alpha_tmp = freq / n_r / const.c_light * np.array([epsilon_x, epsilon_y, epsilon_z, epsilon_ave])
            alpha.append(alpha_tmp)

        alpha = np.array(alpha)

        alpha = alpha / const.Bohrtocm

        f = open(self.path + 'Absorption-' + self.name + '.dat', 'w')  # file name
        f.write(
            "#Erange (eV), Absorption coefficient (cm$^{-1}$) x, y, z, total\n")
        for i in range(len(Erange)):
            f.write("%s %s \n" % (Erange[i], str(alpha[i]).replace('[', ' ').replace(']', ' ')))
        f.close()


    def calculate(self, T: float, density: float):
        prefactor = 4 * self.mu * const.e_charge ** 2 / \
                    (const.hbar ** 2 * const.C ** 3 * const.e_mass ** 2 * self.Vcell_cm)
        occ_e, occ_h = self.get_occ(T, density)
        rates = 1 / float(self.nktot) * occ_e * occ_h * (self.ebc - self.ebv) * (
                self.matrix[:, 0] ** 2 + self.matrix[:, 1] ** 2 + self.matrix[:, 2] ** 2) / 3.0
        rates_x = 1 / float(self.nktot) * occ_e * occ_h * (self.ebc - self.ebv) * (
                self.matrix[:, 0] ** 2 )
        rates_y = 1 / float(self.nktot) * occ_e * occ_h * (self.ebc - self.ebv) * (
                self.matrix[:, 1] ** 2 )
        rates_z = 1 / float(self.nktot) * occ_e * occ_h * (self.ebc - self.ebv) * (
                self.matrix[:, 2] ** 2 )

        rates = prefactor * rates * self.conversion_factor
        rates_x = prefactor * rates_x * self.conversion_factor
        rates_y = prefactor * rates_y * self.conversion_factor
        rates_z = prefactor * rates_z * self.conversion_factor
        return sum(rates), sum(rates_x), sum(rates_y), sum(rates_z)

    def get_fermi_elec_carrier(self, T: float, density_au: float):  # !!!!the fermi coefficients for electrons

        kT = const.kB * T
        df1 = self.Eigs[:, int(min(self.cb)) - 1, :]
        ecbm = df1[:, -1].min()

        efermi1 = -20.0
        efermi2 = 20.0

        density1 = 0.0
        density2 = 0.0

        for ib in range(int(min(self.cb)) - 1, self.nb):
            eig = self.Eigs[:, ib, -1]
            occ1 = 1. / (1 + np.exp((-efermi1 + eig) / kT))
            occ2 = 1. / (1 + np.exp((-efermi2 + eig) / kT))
            density1 = density1 + occ1.sum()
            density2 = density2 + occ2.sum()

        if not self.soc:
            density1 = density1 * 2.0 / self.Vcell_cm
            density2 = density2 * 2.0 / self.Vcell_cm
        else:
            density1 = density1 * 1.0 / self.Vcell_cm
            density2 = density2 * 1.0 / self.Vcell_cm

        if (density1 - density_au) * (density2 - density_au) > 0:
            print("wrong energy window!")
            exit()
        for step in range(100):
            efermi = (efermi1 + efermi2) / 2.0
            density = 0.0

            for ib in range(int(min(self.cb)) - 1, self.nb):
                eig = self.Eigs[:, ib, -1]
                occ = 1. / (1 + np.exp((-efermi + eig) / kT))
                density = density + occ.sum()

            if not self.soc:
                density = density * 2.0 / self.Vcell_cm
            else:
                density = density * 1.0 / self.Vcell_cm

            if abs(density - density_au) / density_au < 1.0E-6:
                break
            if density > density_au:
                efermi2 = efermi
                density2 = density
            else:
                efermi1 = efermi
                density1 = density
            if step >= 99:
                print("bisecting failed!")

        return efermi, ecbm

    def get_fermi_hole_carrier(self, T: float, density_au: float):  # !!!!the fermi coefficients for electrons

        kT = const.kB * T
        df1 = self.Eigs[:, int(max(self.vb)) - 1, :]
        evbm = df1[:, -1].max()

        efermi1 = -20.0
        efermi2 = 20.0

        density1 = 0.0
        density2 = 0.0

        for ib in range(0, int(max(self.vb))):
            eig = self.Eigs[:, ib, -1]
            occ1 = 1. / (1 + np.exp((efermi1 - eig) / kT))
            occ2 = 1. / (1 + np.exp((efermi2 - eig) / kT))
            density1 = density1 + occ1.sum()
            density2 = density2 + occ2.sum()

        if not self.soc:
            density1 = density1 * 2.0 / self.Vcell_cm
            density2 = density2 * 2.0 / self.Vcell_cm
        else:
            density1 = density1 * 1.0 / self.Vcell_cm
            density2 = density2 * 1.0 / self.Vcell_cm

        if (density1 - density_au) * (density2 - density_au) > 0:
            print("wrong energy window!")
            exit()
        for step in range(100):
            efermi = (efermi1 + efermi2) / 2.0
            density = 0.0

            for ib in range(0, int(max(self.vb))):
                eig = self.Eigs[:, ib, -1]
                occ = 1. / (1 + np.exp((efermi - eig) / kT))
                density = density + occ.sum()

            if not self.soc:
                density = density * 2.0 / self.Vcell_cm
            else:
                density = density * 1.0 / self.Vcell_cm

            print(step, efermi, density)

            if abs(density - density_au) / density_au < 1.0E-6:
                break
            if density > density_au:
                efermi1 = efermi
                density1 = density
            else:
                efermi2 = efermi
                density2 = density
            if step >= 99:
                print("bisecting failed!")

        return efermi, evbm

    # calculate quasi fermi level for electrons
    def get_fermi_elec(self, T: float, density_au: float):  # !!!!the fermi coefficients for electrons
        kT = const.kB * T
        efermi1 = -20.0
        efermi2 = 20.0

        eig = np.array(self.ebc)
        ecbm = eig.min()
        kw = 1/self.nktot
        occ1 = 1. / (1 + np.exp((-efermi1 + eig) / kT)) * kw
        occ2 = 1. / (1 + np.exp((-efermi2 + eig) / kT)) * kw
        density1 = occ1.sum() / self.nvb
        density2 = occ2.sum() / self.nvb

        if not self.soc:
            density1 = density1 * 2.0 / self.Vcell_cm
            density2 = density2 * 2.0 / self.Vcell_cm
        else:
            density1 = density1 * 1.0 / self.Vcell_cm
            density2 = density2 * 1.0 / self.Vcell_cm

        if (density1 - density_au) * (density2 - density_au) > 0:
            raise ValueError(" the value of data is out of \
                            the range of density windows")

        for step in range(100):
            efermi = (efermi1 + efermi2) / 2.0
            occ = 1. / (1 + np.exp((-efermi + eig) / kT)) * kw
            density = occ.sum() / self.nvb

            if not self.soc:
                density = density * 2.0 / self.Vcell_cm
            else:
                density = density * 1.0 / self.Vcell_cm

            if abs(density - density_au) / density_au < 1.0E-6:
                break
            if density > density_au:
                efermi2 = efermi
            else:
                efermi1 = efermi
            if step >= 99:
                print("bisecting failed!")

        return efermi, ecbm

    # calculate quasi fermi level for holes
    def get_fermi_hole(self, T, density_au):

        kT = const.kB * T
        efermi1 = -20.0
        efermi2 = 20.0

        eig = self.ebv
        evbm = eig.max()
        kw = 1 / self.nktot
        occ1 = 1. / (1 + np.exp((efermi1 - eig) / kT)) * kw
        occ2 = 1. / (1 + np.exp((efermi2 - eig) / kT)) * kw
        density1 = occ1.sum() / self.ncb
        density2 = occ2.sum() / self.ncb
        if not self.soc:
            density1 = density1 * 2.0 / self.Vcell_cm
            density2 = density2 * 2.0 / self.Vcell_cm
        else:
            density1 = density1 * 1.0 / self.Vcell_cm
            density2 = density2 * 1.0 / self.Vcell_cm

        if (density1 - density_au) * (density2 - density_au) > 0:
            raise ValueError(" the value of data is out of \
                            the range of density windows")

        for step in range(100):
            efermi = (efermi1 + efermi2) / 2.0
            occ = 1. / (1 + np.exp((efermi - eig) / kT)) * kw
            density = occ.sum() / self.ncb
            if not self.soc:
                density = density * 2.0 / self.Vcell_cm
            else:
                density = density * 1.0 / self.Vcell_cm

            if abs(density - density_au) / density_au < 1.0E-6:
                break
            if density > density_au:
                efermi1 = efermi
            else:
                efermi2 = efermi
            if step >= 99:
                print("bisecting failed!")

        return efermi, evbm

    # calculate k coordinates dominate by electrons
    def get_klist_elec(self, T, density_au):
        efermi, ecbm = self.get_fermi_elec_carrier(T, density_au)
        kT = const.kB * T
        for iT in range(100):#, desc="KT"):
            ecutoff = ecbm + iT * kT

            density = 0.0
            klist = np.array([])
            for ib in range(int(min(self.cb))-1, self.nb):
                df1 = self.Eigs[:, ib, :]
                eigs = df1[:, 4]
                kx = df1[:, 1]
                ky = df1[:, 2]
                kz = df1[:, 3]
                for i in range(df1.shape[0]):
                    if eigs[i] < ecutoff:
                        w = 1. / (1 + np.exp((-efermi + eigs[i]) / kT))
                        density = density + w
                        klist = np.append(klist, [ib + 1, kx[i], ky[i], kz[i], w])
            klist = np.resize(klist, (int(len(klist) / 5), 5))
            if not self.soc:
                density = density * 2.0 / self.Vcell_cm
            else:
                density = density * 1.0 / self.Vcell_cm

            print(iT, "kT, density = ", density)
            if abs(density - density_au) / density_au < 1.0E-2:
                print("99% electron density achieved!")
                break

        return efermi, klist, ecutoff

    # calculate k coordinates dominate by holes
    def get_klist_hole(self, T, density_au):
        efermi, evbm = self.get_fermi_hole_carrier(T, density_au)
        kT = const.kB * T
        for iT in range(100):#, desc="KT"):
            ecutoff = evbm - iT * kT
            density = 0.0
            klist = np.array([])
            for ib in range(0, int(max(self.vb))):
                df1 = self.Eigs[:, ib, :]
                eigs = df1[:, 4]
                kx = df1[:, 1]
                ky = df1[:, 2]
                kz = df1[:, 3]
                for i in range(df1.shape[0]):
                    if eigs[i] > ecutoff:
                        w = 1. / (1 + np.exp((efermi - eigs[i]) / kT))
                        density = density + w
                        klist = np.append(klist, [ib + 1, kx[i], ky[i], kz[i], w])
            klist = np.resize(klist, (int(len(klist) / 5), 5))
            if not self.soc:
                density = density * 2.0 / self.Vcell_cm
            else:
                density = density * 1.0 / self.Vcell_cm

            if abs(density - density_au) / density_au < 1.0E-2:
                print("99% electron density achieved!")
                break

        return efermi, klist, ecutoff

    def get_occ(self, T, density_au):
        EF_e, _ = self.get_fermi_elec(T, density_au)
        EF_h, _ = self.get_fermi_hole(T, density_au)

        kT = const.kB * T
        e_e = self.ebc - EF_e
        occ_e = 1. / (1 + np.exp((e_e) / kT))

        e_h = self.ebv - EF_h
        occ_h = 1. / (1 + np.exp((-e_h) / kT))

        return occ_e, occ_h

    def get_Radcoeffi_results(self, T: np.ndarray, density: np.ndarray):
        Rate = []
        for t in tqdm(T, desc="t"):
            for d in density:
                rate, rate_x, rate_y, rate_z = self.calculate(t, d)
                Rate.append(rate/d ** 2)
        Rate = np.array(Rate)
        Bcoeff = Rate.reshape(-1, density.shape[0])
        t1 = 1 / (Bcoeff * density) * 1e9

        f = open(self.path + 'Radiative-' + self.name + '_trail.dat', 'w')  # file name
        f.write(
            "#density (cm^-3),B_%sK (cm^3s^-1), lifetime_%sK (ns)\n" % (T, T))
        for i in range(len(density)):
                f.write("%.5e %s %s \n" % (density[i], str(Bcoeff[:, i]).replace('[',' ').replace(']', ' '),\
                                           str(t1[:, i]).replace('[', ' ').replace(']', ' ')))
        f.close()

        return Bcoeff

    def setup_axes(self, fig, rect):
        tr = Affine2D().skew_deg(30, 0)

        grid_locator = FixedLocator(np.linspace(-0.5, 0.5, 11))
        tick_f = dict()
        for i in np.linspace(-0.5, 0.5, 11):
            tick_f[i] = ''
        for i in np.linspace(-0.5, 0.5, 3):
            tick_f[i] = str(i)
        tick_formatter = DictFormatter(dict(tick_f))

        grid_helper = floating_axes.GridHelperCurveLinear(
            tr, extremes=(-0.5, 0.5, -0.5, 0.5),
            grid_locator1=grid_locator,
            grid_locator2=grid_locator,
            tick_formatter1=tick_formatter,
            tick_formatter2=tick_formatter)

        ax = fig.add_subplot(
            rect, axes_class=floating_axes.FloatingAxes, grid_helper=grid_helper)

        ax.set_xlabel('$k_x (2\pi/a)$')
        ax.set_ylabel('$k_y (2\pi/a)$')

        aux_ax = ax.get_aux_axes(tr)

        return ax, aux_ax

    def plot_k_slice(self, zval: np.ndarray, T: np.ndarray, density: np.ndarray):
        for t in T:
            for d in density:
                for z in zval:
                    self.plot_slice(z, t, d)
                #calculate the k points distribution

    def plot_slice(self, zval: float, T: float, density: float):
        plt.rcParams.update({"font.size": 40})
        efermi_elec, d_k_elec, cutoff_elec = self.get_klist_elec(T, density)
        efermi_hole, d_k_hole, cutoff_hole = self.get_klist_hole(T, density)

        #some correction about the k points coordinates
        k, b, e = self.Eigs.shape
        df1 = self.Eigs.reshape(k*b, e)
        for i in range(1, 4):
            df1[(df1[:, i] > 0.5)][:, i] += -1
        #add the edge of K points
        df = df1
        for i in range(1, 4):
            edge = df[(df[:, i] == -0.5)]
            edge[:, i] = 0.5
            df = np.concatenate((df,edge), axis=0)

        # reorder matrix to make the same band combination together
        order = np.lexsort((df[:, 4], df[:, 3], df[:, 2], df[:, 1]))
        df = df[order].reshape(-1, b, e)

        for i in self.vb:
            # since the band input start by 1 diff from python, so correction is added
            slice_v = df[(df[:, 0, 3] == zval)][:, int(i) - 1, :]
            XX_v = slice_v[:, 1].reshape(self.kmesh[0] + 1, self.kmesh[1] + 1, order='f')
            YY_v = slice_v[:, 2].reshape(self.kmesh[0] + 1, self.kmesh[1] + 1, order='f')
            ZZ0_v = slice_v[:, -1].reshape(self.kmesh[0] + 1, self.kmesh[1] + 1, order='f')

            k_hole = d_k_hole[(d_k_hole[:, 0] == int(i)) & (d_k_hole[:, 3] == zval)]
            kx_hole = k_hole[:, 1]
            ky_hole = k_hole[:, 2]
            w_hole = k_hole[:, -1] # obtain x y coordinates and eigenvalues
            if w_hole.any():
                w_hole = w_hole / min(w_hole) * np.exp(1)

            fig_v = plt.figure(figsize=(20, 15))
            ax_v, aux_ax_v = self.setup_axes(fig_v, 111)
            im_v = aux_ax_v.contourf(XX_v, YY_v, ZZ0_v, 20)
            aux_ax_v.scatter(kx_hole, ky_hole, c='b', s=1 * np.log(w_hole))
            ax_v.set_title('Valence band' + str(int(i)), loc='right')
            fig_v.colorbar(im_v, ax=ax_v)
            plt.savefig('VB' + str(int(i)) +
                        '-' + str(T) + 'K.png')

        for i in self.cb:
            # since the band input start by 1 diff from python, so correction is added
            slice_c = df[(df[:, 0, 3] == zval)][:, int(i) - 1, :]
            XX_c= slice_c[:, 1].reshape(self.kmesh[0] + 1, self.kmesh[1] + 1, order='f')
            YY_c = slice_c[:, 2].reshape(self.kmesh[0] + 1, self.kmesh[1] + 1, order='f')
            ZZ0_c = slice_c[:, -1].reshape(self.kmesh[0] + 1, self.kmesh[1] + 1, order='f')

            k_elec = d_k_elec[(d_k_elec[:, 0] == int(i)) & (d_k_elec[:, 3] == zval)]
            kx_elec = k_elec[:, 1]
            ky_elec = k_elec[:, 2]
            w_elec = k_elec[:, -1]
            if w_elec.any():
                w_elec = w_elec / min(w_elec) * np.exp(1)

            fig_c = plt.figure(figsize=(20, 15))
            ax_c, aux_ax_c = self.setup_axes(fig_c, 111)
            im_c = aux_ax_c.contourf(XX_c, YY_c, ZZ0_c, 20)
            aux_ax_c.scatter(kx_elec, ky_elec, c='r', s=1 * np.log(w_elec))
            ax_c.set_title('Conduction band' + str(int(i)), loc='right')
            fig_c.colorbar(im_c, ax=ax_c)
            plt.savefig('CB' + str(int(i)) + '-' + str(T) + 'K.png')

    def radiative_coeff_plot(self, T: np.ndarray, density: np.ndarray):

        Bcoeff = self.get_Bcoeffi_results(T, density)

        for i in range(T.shape[0]):
            plt.plot(density, Bcoeff[i], 'o-', label='$T$ = %s K' %T[i])

        plt.xscale('log')
        plt.xlim([min(density), max(density)])
        plt.xlabel(r"Carrier density (cm$^{-3}$)")
        plt.ylabel(r"$B$ coefficient (cm$^{3}$s$^{-1}$)")

        plt.legend()
        plt.savefig('Bcoeff-' + self.name + '_trail.png')  # file name
        plt.figure()

        t1 = 1 / (Bcoeff * density) * 1e9
        for i in range(T.shape[0]):
            plt.plot(density, t1[i], 'o-', label='$T$ = %s K' %T[i])

        plt.xscale('log')
        plt.yscale('log')
        plt.xlim([min(density), max(density)])
        plt.xlabel(r"Carrier density (cm$^{-3}$)")
        plt.ylabel(r"Radiative lifetime (ns)")
        plt.legend()
        plt.savefig('lifetime-'+ self.name + '_trail.png')  # file name
        plt.show()



# module for parameter modification
if __name__ == '__main__':
    Density = np.array([1E15, 3.5E15, 1E16, 3.5E16, 1E17, 3.5E17, 1E18, 3.5E18, 1E19])
    Trange = np.array([100, 200, 300])

    T1 = [100,200,300]
    D1 = [1E19]

    mu = 3.9476 #1.48 # refractive index !!!! 1.61 for CsMgBr3 1.58 for CsCaBr3 1.48 for CsMgCl3 3.9476
    Vcell = 183 #394.20 # volume of cell, find in OUTCAR !!! 183.01
    name = 'GaAs'
    Rate = rate_cal(name, mu, Vcell, soc=False, root="../data/GaAs/", start_band=20)

    Emin = 0.5
    Emax = 5.5
    Ne = 301
    Eg = 2.390993
    EgR = 2.31489694402
    VB = [12, 13, 14, 15, 16]#[77, 78, 79, 80]#vbm
    CB = [17, 18, 19, 20, 21]#[81, 82, 83, 84]#cbm

    #Rate.absorption(Emin, Emax, Ne, 0.040)

    #Figure = figure_plot('CsMgBr3', './')
    #Figure.adsorption_coeff_plot(Eg= Eg, EgR= EgR)
    #Figure.radiative_coeff_plot(1E15, 1E19)
    #Figure.plot_k_slice('y', [0.0], T1, D1, VB, CB, True)

    #Rate.plot_k_slice([0.0], T1, D1)
    #Rate.radiative_coeff_plot(Trange, Density)

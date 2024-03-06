# Copyright (c) Xie Zhang
# Distributed under the terms of the MIT License.

import os
import math
import numpy as np
import matplotlib.pyplot as plt
from rate.constants import const
from matplotlib.transforms import Affine2D
from transmatrix_interp.interpolator import bicubic_interp, linear_interp
import mpl_toolkits.axisartist.floating_axes as floating_axes
from mpl_toolkits.axisartist.grid_finder import (FixedLocator, DictFormatter)
from input import read_EIGEN, read_kpt, read_wannier_eigs, read_radiative, read_adsorption

'''Visualization module for all results including carrier distribution

In this module, all the final results under different conditions
can be visualized by output graph. Moreover, the carrier distribution
are also exhibited to help users analyze the mechanism that would 
influence those coefficients.
'''

class figure_plot:
    def __init__(self, name: str, root= './'):
        self.name = name
        self.path = root

        self.Eigs = None
        self.kmesh = None
        self.vb = None
        self.cb = None

    def radiative_coeff_plot(self, xmin, xmax):
        file = self.path + 'Radiative-' + self.name + '_trail.dat'
        T, density, Bcoeff, lifetimes = read_radiative(file)

        for i in range(T.shape[0]):
            plt.plot(density, Bcoeff[i], 'o-', label='$T$ = %s K' % T[i])
        if xmin == None: xmin = min(density)
        if xmax == None: xmax = max(density)
        plt.xscale('log')
        # plt.yscale('log')
        plt.xlim([xmin, xmax])
        #plt.ylim([1e-13, 1e-8])
        plt.xlabel(r"Carrier density (cm$^{-3}$)")
        plt.ylabel(r"$B$ coefficient (cm$^{3}$s$^{-1}$)")

        plt.legend()
        plt.savefig('Bcoeff-' + self.name + '_trail.png')  # file name
        plt.figure()

        t1 = 1 / (Bcoeff * density) * 1e9
        for i in range(T.shape[0]):
            plt.plot(density, t1[i], 'o-', label='$T$ = %s K' % T[i])

        plt.xscale('log')
        plt.yscale('log')
        plt.xlim([xmin, xmax])
        plt.xlabel(r"Carrier density (cm$^{-3}$)")
        plt.ylabel(r"Radiative lifetime (ns)")
        plt.legend()
        plt.savefig('lifetime-' + self.name + '_trail.png')  # file name
        plt.show()

    def adsorption_coeff_plot(self, Eg, Amax = -1, fsize = 12):
        file = self.path + 'Absorption-' + self.name + '.dat'
        data = np.loadtxt(file)
        Erange = data[:, 0]
        alpha = data[:, 1::]
        Emin = Erange.min()
        Emax = Erange.max()
        if Amax == -1 : Amax = 10**math.ceil(math.log(alpha.max(), 10))
        plt.plot(Erange, alpha[:, 3], 'r-', label="average")
        plt.vlines(Eg, 1, Amax, colors='k', linestyles='dashed')
        plt.xlim([Emin, Emax])
        plt.ylim([1, Amax])
        plt.yscale('log')
        plt.text(Eg, 11000, '$E_g$', fontsize=fsize)
        plt.xlabel("Energy (eV)", fontsize=fsize)
        plt.ylabel("Absorption coefficient (cm$^{-1}$)", fontsize=fsize)
        plt.xticks(fontsize=fsize)
        plt.yticks(fontsize=fsize)
        plt.savefig(self.path + "Absorption-" + self.name + ".pdf")

    def plot_k_slice(self, direct: str, zval: np.ndarray, Vcell: float, T: np.ndarray, density: np.ndarray,\
                     VB: np.ndarray, CB: np.ndarray, soc=False, startband = 0, fontsize = 38):
        self.Vcell_cm = Vcell * const.Ang2cm ** 3
        FineKptfile = os.path.join(self.path + 'wannier90_geninterp.kpt')
        FineEigfile = os.path.join(self.path + 'wannier90_geninterp.dat')
        Ori_eigfile = os.path.join(self.path + 'EIGENVAL')
        self.soc = soc
        nktot, kpts = read_kpt(FineKptfile)
        ksize, _, self.vborder, Eigenvalue = read_EIGEN(Ori_eigfile, self.soc)
        self.interpolator = linear_interp(1)
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

        for i in range(kpts.shape[0]):
            self.Eigs[i, :, 1:4] = kpts[i, 1:4]

        Eigs = self.Eigs[:, :, -1]
        evbm = max(Eigs[:, self.vborder-1 - startband])
        self.Eigs[:, :, -1] -= evbm
        self.vb = [vb - startband for vb in VB]
        self.cb = [cb - startband for cb in CB]
        self.nb = int(eignum / nktot)

        for t in T:
            for d in density:
                for z in zval:
                    print(self.kmesh)
                    self.plot_slice(direct, z, t, d, fontsize, startband)
                #calculate the k points distribution

    def plot_slice(self, direct, zval: float, T: float, density: float, fontsize:float, startband:int):
        plt.rcParams.update({"font.size": fontsize})
        efermi_elec, d_k_elec, cutoff_elec = self.get_klist_elec(T, density)
        efermi_hole, d_k_hole, cutoff_hole = self.get_klist_hole(T, density)

        #some correction about the k points coordinates
        k, b, e = self.Eigs.shape
        df1 = self.Eigs.reshape(k*b, e)
        for i in range(1, 4):
            df1[(df1[:, i] > 0.5)][:, i] += -1
        #add the edge of K points
        df = df1
        kmesh = self.kmesh.copy()
        for i in range(1, 4):
            edge = df[(df[:, i] == 0.5)]
            edge[:, i] = -0.5
            if edge.shape[0] > 0:
                kmesh[i-1] = kmesh[i-1] + 1
            df = np.concatenate((df, edge), axis=0)

        # reorder matrix to make the same band combination together
        order = np.lexsort((df[:, 4], df[:, 3], df[:, 2], df[:, 1]))
        df = df[order].reshape(-1, b, e)

        coor = ['x', 'y', 'z']
        coor.remove(direct)
        cmap = {'x': 1, 'y': 2, 'z': 3} # map the direction string

        for i in self.vb:
            # since the band input start by 1 diff from python, so correction is added
            slice_v = df[(df[:, 0, cmap[direct]] == zval)][:, int(i) - 1, :]
            XX_v = slice_v[:, cmap[coor[0]]].reshape(kmesh[cmap[coor[0]]-1], kmesh[cmap[coor[1]]-1], order='C')
            YY_v = slice_v[:, cmap[coor[1]]].reshape(kmesh[cmap[coor[0]]-1], kmesh[cmap[coor[1]]-1], order='C')
            ZZ0_v = slice_v[:, -1].reshape(kmesh[cmap[coor[0]]-1], kmesh[cmap[coor[1]]-1], order='C')

            k_hole = d_k_hole[(d_k_hole[:, 0] == int(i)) & (d_k_hole[:, cmap[direct]] == zval)]
            kx_hole = k_hole[:, cmap[coor[0]]]
            ky_hole = k_hole[:, cmap[coor[1]]]
            w_hole = k_hole[:, -1] # obtain x y coordinates and eigenvalues
            if w_hole.any():
                w_hole = w_hole / min(w_hole) * np.exp(1)

            fig_v = plt.figure(figsize=(30, 20), dpi=100)
            ax_v, aux_ax_v = self.setup_axes(coor, fig_v, 121)
            #im_v = aux_ax_v.contourf(XX_v, YY_v, ZZ0_v, 30, cmap='coolwarm', alpha=0.6)
            #aux_ax_v.scatter(kx_hole, ky_hole, c='yellow', edgecolors='salmon', linewidths=0, s=42, alpha=0.42)
            im_v = aux_ax_v.contourf(XX_v, YY_v, ZZ0_v, 20)
            aux_ax_v.scatter(kx_hole, ky_hole, c='b', s=10 * np.log(w_hole))
            ax_v.set_title('VB' + str(int(abs(i-np.max(self.vb))+1)), loc='right')
            fig_v.colorbar(im_v, ax=ax_v, fraction=0.032, ticks=[np.min(ZZ0_v), 0.5*(np.min(ZZ0_v)+np.max(ZZ0_v)), np.max(ZZ0_v)])
            plt.savefig(self.path + self.name + '_' + direct + '=' + str(zval) + \
                        '-VB' + str(int(i+startband)) + '-' + str(density) + '-' + str(T) + 'K.png')

        for i in self.cb:
            # since the band input start by 1 diff from python, so correction is added
            slice_c = df[(df[:, 0, cmap[direct]] == zval)][:, int(i) - 1, :]
            XX_c= slice_c[:, cmap[coor[0]]].reshape(kmesh[cmap[coor[0]]-1], kmesh[cmap[coor[1]]-1], order='C')
            YY_c = slice_c[:, cmap[coor[1]]].reshape(kmesh[cmap[coor[0]]-1], kmesh[cmap[coor[1]]-1], order='C')
            ZZ0_c = slice_c[:, -1].reshape(kmesh[cmap[coor[0]]-1], kmesh[cmap[coor[1]]-1], order='C')

            k_elec = d_k_elec[(d_k_elec[:, 0] == int(i)) & (d_k_elec[:, cmap[direct]] == zval)]
            kx_elec = k_elec[:, cmap[coor[0]]]
            ky_elec = k_elec[:, cmap[coor[1]]]
            w_elec = k_elec[:, -1]
            if w_elec.any():
                w_elec = w_elec / min(w_elec) * np.exp(1)

            fig_c = plt.figure(figsize=(30, 20), dpi=100)
            ax_c, aux_ax_c = self.setup_axes(coor, fig_c, 121)
            im_c = aux_ax_c.contourf(XX_c, YY_c, ZZ0_c, 20)
            aux_ax_c.scatter(kx_elec, ky_elec, c='r', s=10 * np.log(w_elec))
            ax_c.set_title('CB' + str(int(i-np.min(self.cb)+1)), loc='right')
            fig_c.colorbar(im_c, ax=ax_c, fraction=0.032, ticks=[np.min(ZZ0_c), 0.5*(np.min(ZZ0_c)+np.max(ZZ0_c)), np.max(ZZ0_c)])
            plt.savefig(self.path + self.name + '_' +direct + '=' + str(zval) + \
                        '-CB' + str(int(i+startband)) + '-' + str(density) + '-' + str(T) + 'K.png')

    # calculate quasi fermi level for holes
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

    #calculate quasi fermi level for electrons
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

    # calculate k coordinates dominate by electrons
    def get_klist_elec(self, T, density_au):
        efermi, ecbm = self.get_fermi_elec_carrier(T, density_au)
        kT = 8.6173303E-5 * T
        for iT in range(100):#, desc="KT"):
            ecutoff = ecbm + iT * kT

            density = 0.0
            klist = np.empty((0, 5))
            for ib in range(int(min(self.cb))-1, self.nb):
                df1 = self.Eigs[:, ib, :]
                eigs = df1[:, 4]
                kx = df1[:, 1]
                ky = df1[:, 2]
                kz = df1[:, 3]
                order = eigs < ecutoff
                w = 1. / (1 + np.exp((-efermi + eigs[order]) / kT))
                density = density + sum(w)
                first_col = np.ones(kx[order].shape[0]) * (ib + 1)
                klist_tmp = np.column_stack((first_col, kx[order], ky[order], kz[order], w))
                klist = np.concatenate((klist, klist_tmp), axis=0)

            #klist = np.resize(klist, (int(len(klist) / 5), 5))
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
        kT = 8.6173303E-5 * T
        for iT in range(100):#, desc="KT"):
            ecutoff = evbm - iT * kT
            density = 0.0
            klist = np.empty((0, 5))
            for ib in range(0, int(max(self.vb))):
                df1 = self.Eigs[:, ib, :]
                eigs = df1[:, 4]
                kx = df1[:, 1]
                ky = df1[:, 2]
                kz = df1[:, 3]
                order = eigs > ecutoff
                w = 1. / (1 + np.exp((efermi - eigs[order]) / kT))
                density = density + sum(w)
                first_col = np.ones(kx[order].shape[0]) * (ib + 1)
                klist_tmp = np.column_stack((first_col, kx[order], ky[order], kz[order], w))
                klist = np.concatenate((klist, klist_tmp), axis=0)

            #klist = np.resize(klist, (int(len(klist) / 5), 5))
            if not self.soc:
                density = density * 2.0 / self.Vcell_cm
            else:
                density = density * 1.0 / self.Vcell_cm

            print(iT, "kT, density = ", density)
            if abs(density - density_au) / density_au < 1.0E-2:
                print("99% electron density achieved!")
                break

        return efermi, klist, ecutoff

    #parameter settings for axes in graph
    def setup_axes(self, coor, fig, rect):
        tr = Affine2D().skew_deg(0, 0)

        grid_locator = FixedLocator(np.around(np.linspace(-0.40, 0.40, 5), 2))
        tick_f = dict()
        for i in np.linspace(-0.40, 0.40, 5):
            i = round(i, 2)
            tick_f[i] = ''
        for i in np.linspace(-0.40, 0.40, 5):
            i = round(i, 2)
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

        ax.set_xlabel('$k_' + coor[0] + '(2\pi/a)$')
        ax.set_ylabel('$k_' + coor[1] + '(2\pi/a)$')

        aux_ax = ax.get_aux_axes(tr)

        return ax, aux_ax

# module for parameter modification
if __name__ == '__main__':
    T1 = [100, 300]
    D1 = [1E19]
    VB = [34, 35, 36]#[77, 78, 79, 80]#vbm
    CB = [37, 38, 39]#[81, 82, 83, 84]#cbm

    Eg = 8.0
    EgR = 1.424

    Vcell = 180
    Figure = figure_plot('GaAs', '../test_files/')
    #Figure.adsorption_coeff_plot(Eg=Eg, EgR=EgR)
    #Figure.radiative_coeff_plot(1E15, 1E19)
    Figure.plot_k_slice('z', [0.0], Vcell, T1, D1, VB, CB, False, 20)
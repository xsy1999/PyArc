# Copyright (c) Xie Zhang
# Distributed under the terms of the MIT License.

import os
import numpy as np
from transmatrix_interp.interpolator import bicubic_interp, linear_interp
from tqdm import tqdm
from input import read_transmatrix, read_kpt, read_EIGEN

'''The Module for eigenvalue and matrix value interpolation

In this module both eigenvalues and matrix values are interpolated
through the interpolation scheme in interpolator module. And
the results are output as files for check and later calculation
'''

class Transmatrix:
    """parent class for all interpolators

    some initial process and file writing process are achieved here

    interpolate method could be changed according to different method
    """
    def __init__(self, Nvbn: np.ndarray, Ncbn: np.ndarray, root='./', method = ['linear', 'cubic']):
        self.path = root
        self.ncb = Ncbn
        self.nvb = Nvbn
        IBZKPT = os.path.join(self.path + 'IBZKPT')
        Transpath = os.path.join(self.path + 'Transmatrix')
        matrixpath = os.path.join(self.path + 'matrix.npy')
        Ori_eigfile = os.path.join(self.path + 'EIGENVAL')
        knum, kpts = read_kpt(IBZKPT)
        self.kpts = kpts[:, 0:-1] # The kpoints index of original mesh
        self.knum = knum #number of k points
        np.save("kpts.npy", kpts)
        if os.path.exists(matrixpath):
            os.remove(matrixpath)

        ksize, _, self.vborder, self.Eigenvalue = read_EIGEN(Ori_eigfile, False)

        matrix = []
        if os.path.exists(Transpath):
            Transnum, Trans = read_transmatrix(Transpath)
        else: raise IOError("No trans matrix file")

        if not Transnum % knum == 0:  # check if whether the trans file match the IBZ file
            raise ValueError("Error, the Tranmatrix file does not match the IBZKPT file")

        # calculate the trans values of xyz direction and insert their kpoints coordinates
        for i in tqdm(range(Transnum), desc='Reoganizing'):
            if int(Trans[i][1]) in Ncbn and int(Trans[i][2]) in Nvbn:
                kindex = int(i // (Transnum / knum))
                trans = [int(Trans[i][1]), int(Trans[i][2]), float(kpts[kindex, 0]),\
                         float(kpts[kindex, 1]), float(kpts[kindex, 2]),\
                         (float(Trans[i][5]) ** 2 + float(Trans[i][6]) ** 2) ** 0.5,\
                         (float(Trans[i][7]) ** 2 + float(Trans[i][8]) ** 2) ** 0.5,\
                         (float(Trans[i][9]) ** 2 + float(Trans[i][10]) ** 2) ** 0.5]
                matrix.append(trans)
        matrix = np.array(matrix)
        order = np.lexsort((matrix[:, 0], matrix[:, 1]))  # reorder matrix to make the same band combination together
        matrix = matrix[order]
        self.matrix = matrix
        self.bc = int(matrix.shape[0] // kpts.shape[0])
        self.interpolator = method
        np.save("matrix.npy", matrix)

    def write_to_file(self, magx = 1, magy = 1, magz = 1) -> None:
        if magx*magy*magz < 1:
            raise ValueError("the magnification value is wrong")
        nvb = self.nvb.shape[0]
        matrix = self.matrix
        bc = self.bc  # number of different band combination
        cbs = self.ncb
        vbs = self.nvb
        order = np.lexsort((matrix[:, 4], matrix[:, 3], matrix[:, 2], matrix[:, 1], matrix[:, 0]))
        matrix = matrix[order]

        if self.interpolator[0] == 'linear':
            interpolator = linear_interp(1)
        elif self.interpolator[0] == 'cubic':
            interpolator = bicubic_interp(1)

        # ensure the size of original reciprocal k mesh
        ncx = np.unique(np.around(matrix[:, 2], 4)).shape[0]
        ncy = np.unique(np.around(matrix[:, 3], 4)).shape[0]
        ncz = np.unique(np.around(matrix[:, 4], 4)).shape[0]

        mkptx = max(matrix[:, 2])
        mkpty = max(matrix[:, 3])
        mkptz = max(matrix[:, 4])

        nintx = int(magx * ncx)
        ninty = int(magy * ncy)
        nintz = int(magz * ncz)

        mat_mesh = matrix[:, 5:8].reshape(-1, ncx, ncy, ncz, 3)

        # generate fine k mesh
        ff_1dx = [i + 1.0 / nintx for i in np.linspace(mkptx - 1, mkptx, nintx, endpoint=False)]
        ff_1dy = [i + 1.0 / ninty for i in np.linspace(mkpty - 1, mkpty, ninty, endpoint=False)]
        ff_1dz = [i + 1.0 / nintz for i in np.linspace(mkptz - 1, mkptz, nintz, endpoint=False)]
        XXf, YYf, ZZf = np.meshgrid(ff_1dx, ff_1dy, ff_1dz, indexing='ij')

        ############################################################
        #output dense mesh trans matrix value
        matpath = self.path + 'matrix_fine.dat'
        f = open(matpath, 'w')
        # for each band combination
        for i in tqdm(range(bc), desc='Bands Combination:'):
            mat_f = interpolator.interpolate_kmesh(mat_mesh[i, :] ** 2, magx, magy, magz)
            #matrix value after interpolation

            cb = int(cbs[int(i // nvb)]) # relavent cb
            vb = int(vbs[int(i % nvb)]) # relavent vb

            mat_f = np.sqrt(mat_f)
            for kx, ky, kz, mx, my, mz in \
                    zip(XXf.flatten(), YYf.flatten(), ZZf.flatten(), mat_f[:,0], mat_f[:,1], mat_f[:,2]):#[0], matx_f[1], matx_f[2]):
                #ensure the range of k coordinates in three directions
                kx = 0.5 - ((0.5 - kx) % 1)
                if abs(kx + 0.5) < 1e-5:
                    kx += 1
                ky = 0.5 - ((0.5 - ky) % 1)
                if abs(ky + 0.5) < 1e-5:
                    ky += 1
                kz = 0.5 - ((0.5 - kz) % 1)
                if abs(kz + 0.5) < 1e-5:
                    kz += 1
                f.write("%s %s %.4f %.4f %.4f %.4f %.4f %.4f\n" \
                        % (cb, vb, kx, ky, kz, mx, my, mz))
        f.close()

        ############################################################
        #output kpoints coordinates for wannier interpolation
        kptpath = self.path + 'wannier90_geninterp.kpt'
        f = open(kptpath, 'w')
        f.write('kmesh %s %s %s\n' % (nintx, ninty, nintz))
        f.write('crystal\n')
        f.write('%s\n' % int(nintx * ninty * nintz))
        kidx = 1 # num of k points calculated
        for kx, ky, kz in zip(XXf.flatten(), YYf.flatten(), ZZf.flatten()):
            kx = 0.5 - ((0.5 - kx) % 1)
            if abs(kx + 0.5) < 1e-5:
                kx += 1
            ky = 0.5 - ((0.5 - ky) % 1)
            if abs(ky + 0.5) < 1e-5:
                ky += 1
            kz = 0.5 - ((0.5 - kz) % 1)
            if abs(kz + 0.5) < 1e-5:
                kz += 1
            f.write("%s %.3f %.3f %.3f\n" % (kidx, kx, ky, kz))
            kidx += 1
        f.close()
        print("kpoints for interpolating are generated, please use eigenvalue interpolators to generate fine eigenvalue mesh")

        ############################################################
        # output eigenvalues file by default cubic or linear method
        if self.interpolator[1] == 'linear':
            eigen_interpolator = linear_interp(1)
        elif self.interpolator[1] == 'cubic':
            eigen_interpolator = bicubic_interp(1)
        #EIGEN = eigen_interpolator.interpolate_kmesh(self.Eigenvalue, magx, magy, magz)
        Eigpath = os.path.join(self.path + 'Eigen_geninterp.dat')
        f = open(Eigpath, 'w')
        f.write('# Written data \n')
        f.write('# Input file comment: kmesh %s %s %s\n' % (nintx, ninty, nintz))
        f.write('#  Kpt_idx   K_x (1/ang)   K_y (1/ang)   K_z (1/ang)   Energy (eV)\n')
        kidx = 1 #num of k points calculated
        for kx, ky, kz in zip(XXf.flatten(), YYf.flatten(), ZZf.flatten()):
            kx = 0.5 - ((0.5 - kx) % 1)
            if abs(kx + 0.5) < 1e-5:
                kx += 1
            ky = 0.5 - ((0.5 - ky) % 1)
            if abs(ky + 0.5) < 1e-5:
                ky += 1
            kz = 0.5 - ((0.5 - kz) % 1)
            if abs(kz + 0.5) < 1e-5:
                kz += 1
            eig = eigen_interpolator.interpolate_kpoints(self.Eigenvalue, np.array([kx,ky,kz])).squeeze()
            for bnd in range(eig.shape[0]):
                f.write("%s   %.3f   %.3f   %.3f   %.5f\n" % (kidx, kx, ky, kz, eig[bnd]))
            kidx += 1
        f.close()

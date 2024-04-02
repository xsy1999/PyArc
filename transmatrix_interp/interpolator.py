# Copyright (c) Xie Zhang
# Distributed under the terms of the MIT License.

import math

import numpy as np
from scipy.interpolate import RegularGridInterpolator

'''the interpolation scheme module

In this module, linear and cubic interpolation have
been achieved for interpolating k sampling mesh
'''
############################################################

class base_interpolator:
	def __init__(self, m=1):
		self.m = 1
	def interpolate(self, nf_x: int, nf_y: int, nf_z: int):
		raise NotImplementedError

#adjustable part
class linear_interp(base_interpolator):
	def __init__(self, m) -> None:
		super().__init__(m)

	def interpolate_kmesh(self, data: np.ndarray, magx: int, magy:  int, magz: int):
		'''
		if self.dimension + [nf_x, nf_y, nf_z].count(1) != 3: #check dimension
			raise ValueError("enlarge parameters are not compatable\
			with dimension of the input mesh")
		'''
		nc_x = data.shape[0]
		nc_y = data.shape[1]
		nc_z = data.shape[2]

		nf_x = int(nc_x * magx)
		nf_y = int(nc_y * magy)
		nf_z = int(nc_z * magz)

		dataf = []

		for i in range(nf_x):
			nx1 = int(((i + 1) * nc_x) // nf_x) % nc_x
			nx0 = nx1 - 1
			xd = (((i + 1) * nc_x) % nf_x) / nf_x

			for j in range(nf_y):
				ny1 = int(((j + 1) * nc_y) // nf_y) % nc_y
				ny0 = ny1 - 1
				yd = (((j + 1) * nc_y) % nf_y) / nf_y

				for k in range(nf_z):
					nz1 = int(((k + 1) * nc_z) // nf_z) % nc_z
					nz0 = nz1 - 1
					zd = (((k + 1) * nc_z) % nf_z) / nf_z

					c000 = data[nx0, ny0, nz0] * (1 - xd) * (1 - yd) * (1 - zd)
					c001 = data[nx0, ny0, nz1] * (1 - xd) * (1 - yd) * zd
					c010 = data[nx0, ny1, nz0] * (1 - xd) * yd * (1 - zd)
					c011 = data[nx0, ny1, nz1] * (1 - xd) * yd * zd
					c100 = data[nx1, ny0, nz0] * xd * (1 - yd) * (1 - zd)
					c101 = data[nx1, ny0, nz1] * xd * (1 - yd) * zd
					c110 = data[nx1, ny1, nz0] * xd * yd * (1 - zd)
					c111 = data[nx1, ny1, nz1] * xd * yd * zd

					c = c000 + c001 + c010 + c011 + c100 + c101 + c110 + c111
					dataf.append(c)
		return dataf

	def interpolate_kpoints(self, data: np.ndarray, kpoints: np.ndarray):
		nc_x = data.shape[0]
		nc_y = data.shape[1]
		nc_z = data.shape[2]

		meshx = np.array([i / nc_x for i in range(nc_x)]) - (math.ceil(nc_x / 2) - 1)/nc_x
		meshy = np.array([i / nc_y for i in range(nc_y)]) - (math.ceil(nc_y / 2) - 1) / nc_y
		meshz = np.array([i / nc_z for i in range(nc_z)]) - (math.ceil(nc_z / 2) - 1) / nc_z

		dataf = []

		for kpoint in kpoints:
			nx0 = np.sum(meshx < kpoint[0]) - 1
			nx1 = nx0 + 1
			xd = min(abs(kpoint[0] - meshx[nx0]), abs(kpoint[0] - meshx[nx0] + 1))

			ny0 = np.sum(meshy < kpoint[1]) - 1
			ny1 = ny0 + 1
			yd = min(abs(kpoint[0] - meshx[nx0]), abs(kpoint[0] - meshx[nx0] + 1))

			nz0 = np.sum(meshz < kpoint[2]) - 1
			nz1 = nz0 + 1
			zd = min(abs(kpoint[0] - meshx[nx0]), abs(kpoint[0] - meshx[nx0] + 1))

			c000 = data[nx0, ny0, nz0] * (1 - xd) * (1 - yd) * (1 - zd)
			c001 = data[nx0, ny0, nz1] * (1 - xd) * (1 - yd) * zd
			c010 = data[nx0, ny1, nz0] * (1 - xd) * yd * (1 - zd)
			c011 = data[nx0, ny1, nz1] * (1 - xd) * yd * zd
			c100 = data[nx1, ny0, nz0] * xd * (1 - yd) * (1 - zd)
			c101 = data[nx1, ny0, nz1] * xd * (1 - yd) * zd
			c110 = data[nx1, ny1, nz0] * xd * yd * (1 - zd)
			c111 = data[nx1, ny1, nz1] * xd * yd * zd

			c = c000 + c001 + c010 + c011 + c100 + c101 + c110 + c111
			dataf.append(c)

		return np.array(dataf)

############################################################

class bicubic_interp(base_interpolator):
	def __init__(self, m) -> None:
		super().__init__(m)

	def interpolate_kmesh(self, data: np.ndarray, magx: int, magy: int, magz: int):
		'''
		if self.dimension + [nf_x, nf_y, nf_z].count(1) != 3: #check dimension
			raise ValueError("enlarge parameters are not compatable\
			with dimension of the input mesh")
		'''
		nc_x = data.shape[0]
		nc_y = data.shape[1]
		nc_z = data.shape[2]

		nf_x = int(data.shape[0] * magx)
		nf_y = int(data.shape[1] * magy)
		nf_z = int(data.shape[2] * magz)

		data = np.pad(data, ((1,1), (1,1), (1,1), (0,0)), 'wrap')

		sr = []
		x = np.arange(-1, nc_x+1, 1)
		y = np.arange(-1, nc_y+1, 1)
		z = np.arange(-1, nc_z+1, 1)
		x_new = np.arange(0, nc_x, 1/magx) if nf_x % 2 == 0 else np.arange(0, nc_x, 1/magx) - 1/magx
		y_new = np.arange(0, nc_y, 1/magy) if nf_y % 2 == 0 else np.arange(0, nc_y, 1/magy) - 1/magy
		z_new = np.arange(0, nc_z, 1/magz) if nf_z % 2 == 0 else np.arange(0, nc_z, 1/magz) - 1/magz
		xg, yg, zg = np.meshgrid(x_new, y_new, z_new)
		f = RegularGridInterpolator((x, y, z), data, method='cubic')
		for x1, y1, z1 in zip(xg.flatten(), yg.flatten(), zg.flatten()):
			result = np.squeeze(f([x1,y1,z1]))
			sr.append(result)
		return np.array(sr)

	def interpolate_kpoints(self, data: np.ndarray, kpoints: np.ndarray):
		nc_x = data.shape[0]
		nc_y = data.shape[1]
		nc_z = data.shape[2]
		data = np.pad(data, ((1,1), (1,1), (1,1), (0,0)), 'wrap')

		x = (np.arange(0, nc_x + 2, 1) - math.ceil(nc_x/2))/nc_x
		y = (np.arange(0, nc_y + 2, 1) - math.ceil(nc_y/2))/nc_y
		z = (np.arange(0, nc_z + 2, 1) - math.ceil(nc_z/2))/nc_z
		f = RegularGridInterpolator((x, y, z), data, method='cubic')


		sr = f(kpoints)
		return sr

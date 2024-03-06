# Copyright (c) Xie Zhang
# Distributed under the terms of the MIT License.

'''Constant module

In the constant module, all the physical parameters like light speed
are stored and prepared for subsequent calculations.
'''

class Const(object):
    class ConstError(TypeError):
        pass

    class ConstCaseError(ConstError):
        pass

    def __setattr__(self, name, value):
        if name in self.__dict__: # the name of an attribute should not exist
            raise self.ConstError("Can't change const.%s" % name)

        self.__dict__[name] = value

const = Const()
const.e_charge = 4.80332E-10
e_charge = 4.80332E-10

const.hbar = 1.0545716E-27
const.e_mass = 9.10938215E-28
const.abohr = 5.2917720859E-9
const.C = 2.99792458E10
const.eV2erg = 1.602E-12
const.Ang2cm = 1E-8
const.ha2eV = 27.2114
const.kB = 8.6173303E-5
const.eVtoHar = 1.0 / 27.2114 # eV to Hartree conversion
const.AngtoBohr = 1.0 / 0.52917720859 # Angstrom to Bohr radius
const.c_light = 137.035999070 # speed of light in Hartree units
const.Bohrtocm = 1 / const.AngtoBohr * 1E-8
const.sigma = 0.012 * const.eVtoHar
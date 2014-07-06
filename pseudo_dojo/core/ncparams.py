from __future__ import print_function, division

from collections import namedtuple


class NcProjector(namedtuple("NCPROJ", "n l rcut scheme")):
    """
    Descriptor for norm-conserving projectors.

    .. attributes:

        n:
            Principal quantum number (the one associated to the 
            AE wavefunction that has been pseudized.
        l:
            Angular momentum.
        rcut:
            Cutoff radius in Bohr.
        scheme:
            String defining the pseudization scheme. 
    """
    @classmethod
    def asprojector(cls, obj):
        if isinstance(obj, cls):
            return obj
        else:
            return cls(**obj)

#class NLCC(dict):

class NcParams(object):
    """
    This object gathers the parameters used to generate a norm-conserving pseudo.

    .. attributes:

        reference_conf:
            String defining the reference configuration.
        Z:
            Nuclear Charge.
        Z_val:
            Number of valence electrons.
        l_max:
            Maximum angular momentum used in the pseudization procedures.
        l_local:
            Angular momemtum used for the local part.
        projectors:
            List of `NcProjector` instances.
        nlcc:
            Defines the treatment of the non-linear core-correction.
        wave_equation:
            String defining the type of Hamiltonian. Possible values:
            (schrodinger, scalar-relativistic, relativistic).
        xc_functional
            String defining the exchange-correlation functional (APE notation).
    """
    def __init__(self, reference_conf, Z, Z_val, l_max, l_local, projectors, 
                 nlcc, wave_equation=None, xc_functional=None, **extra_kwargs):

        self.reference_conf = reference_conf
        self.Z = Z
        self.Z_val = Z_val
        self.l_max = l_max
        self.l_local = l_local
        self.projectors = projectors
        self.nlcc = nlcc
        self.wave_equation = wave_equation
        self.xc_functional = xc_functional
        self._extra_kwargs = extra_kwargs

    @property
    def has_nlcc(self):
        return bool(self.nlcc)

    def set_wave_equation(self, wave_equation):
        self.wave_equation = wave_equation

    def set_xc_functional(self, xc_functional):
        self.xc_functional = xc_functional

    #@classmethod
    #def from_dict(cls, d):

    #@property
    #def to_dict(self):

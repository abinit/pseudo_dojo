"""
systemcodes are (ab inito) codes that can be used to calculate a specific property of a real system, crystal or
molecule. They can be either all electron or pp
"""


__author__ = 'setten'

from abc import ABCMeta, abstractproperty, abstractmethod


class AbstractSystemCodeInterface(object):
    """
    Abstract class for system code interfaces.
    here we should list all methods that are used in the tests
    """
    __metaclass__ = ABCMeta

    @abstractproperty
    def is_all_electron(self):
        """
        return True if the code is all electron, False if it is pseudo
        """

    @abstractproperty
    def periodicity(self):
        """
        return a tuple of booleans, that are true is the code can do calculations in this dimension
        """

    @abstractmethod
    def write_pseudo(self, path, data_set):
        """
        method to write a pp file for this code to path defined by data_set
        """

    @abstractmethod
    def calculate_total_energy(self, structure, kpoints, energy_cutoff, test):
        """
        method to calculate and return the total energy of the defined system and parameters
        if test is True, no calculation should be performed only True should be returned if this method is available
        for this code
        """


class Abinit(AbstractSystemCodeInterface):
    """
    concrete implementation for abinit
    """
    def __init__(self):
        self._periodicity = (False, False, False, True)
        self._is_all_electron = False

    @property
    def is_all_electron(self):
        return self._is_all_electron

    @property
    def periodicity(self):
        return self._periodicity

    def calculate_total_energy(self, structure, kpoints, energy_cutoff, test):
        if test:
            return False
        raise NotImplementedError

    def write_pseudo(self, path, data_set):
        raise NotImplementedError


def get_system_code_interface(code):
    """
    factory function to return an instance of a systems code interface
    """
    cls = {'abinit': Abinit}
    return cls[code]()
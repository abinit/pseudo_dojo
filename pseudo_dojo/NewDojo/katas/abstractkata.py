__author__ = 'setten'

from abc import ABCMeta, abstractproperty, abstractmethod
from pseudo_dojo.NewDojo.codes.systemcodes import get_system_code_interface


class AbstractKata(object):
    """
    Abstract class for defining the interface for kata's (tests)
    the dojo should only use this interface in implementing
    """

    def set_code_interface(self, code):
        self.code_interface = get_system_code_interface(code)


class DeltaFactor(AbstractKata):
    """

    """


def get_kata(kata):
    """
    factory function to return an instance of a systems code interface
    """
    cls = {'deltafactor': DeltaFactor()}
    return cls[kata]()
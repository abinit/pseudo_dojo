"""
generator coded are coded that can be used to create a pp from an input set
"""


__author__ = 'setten'

from abc import ABCMeta, abstractproperty, abstractmethod


class AbstractGeneratorInterface(object):
    """
    Here we define all the methods that all generator interfaces should provide.
    We may distiguish here between really mandatory and not so mandatory. the lather should still be implemented
    by concrete implementations but returning a None we called
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def create(self, definition):
        """
        This method should take the definition, see classes in pseudos, and return a data_set, see classes in pseudos
        """


class ApeGeneratorInterface(AbstractGeneratorInterface):
    """
    The interface for APE as a pp generator
    """

    def create(self, definition):
        raise NotImplementedError


######
# API
######


def get_generator_interface(code):
    """
    factory function to return a instance of a generator code interface
    """
    cls = {'ape': ApeGeneratorInterface}
    return cls[code]()
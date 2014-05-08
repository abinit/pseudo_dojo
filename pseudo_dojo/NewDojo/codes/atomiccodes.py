"""
atomic codes are codes that can be used to do single atom calculations on grid, both all electron and
using a pp
"""

__author__ = 'setten'

from abc import ABCMeta, abstractproperty, abstractmethod


class AbstractAtomicCodeIntrerface(object):
    """
    """
    __metaclass__ = ABCMeta


class ApeAtomicInterface(AbstractAtomicCodeIntrerface):
    """
    """


######
# API
######


def get_atomic_code_interface(code):
    """
    factory function to return a instance of an atomic code interface
    """
    cls = {'ape': ApeAtomicInterface}
    return cls[code]()
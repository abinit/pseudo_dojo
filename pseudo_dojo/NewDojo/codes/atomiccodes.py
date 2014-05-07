"""
atomic codes are codes that can be used to do single atom calculations on grid, both all electron and
using a pp
"""

__author__ = 'setten'

from abc import ABCMeta, abstractproperty, abstractmethod

def get_system_code_interface(code):
    """
    factory function to return a instance of an atomic code interface
    """
    cls = {}
    return cls[code]()
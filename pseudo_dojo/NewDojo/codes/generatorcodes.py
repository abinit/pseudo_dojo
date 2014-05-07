"""
generator coded are coded that can be used to create a pp from an input set
"""


__author__ = 'setten'

from abc import ABCMeta, abstractproperty, abstractmethod


def get_generator_interface(code):
    """
    factory function to return a instance of a generator code interface
    """
    cls = {}
    return cls[code]()
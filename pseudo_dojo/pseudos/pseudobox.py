__author__ = 'setten'

import abc
import collections


class Table(collections.set):
    """
    list of pseudo boxes
    """
    def __init__(self, input):
        if input is list of boxes
            pass
        elif input is string
            make

    @property
    def version(self):
        """
        version of the table
        """

    @property
    def xc_type(self):
        """
        readable string e.g. PBE, ...
        """

    @property
    def pp_type(self):
        """
        readable string e.g PAW, NC, USPP,
        """

    @property
    def pp_generator(self):
        """
        code used to generate the PP, and version
        """

    def consistency_check(self):
        """
        test if the elements comply to the properties of the table
        """


class AbstractPseudoBox(object):
    """
    concrte implentations for pseudo generators
    """
    def __init__(self):


    @abc.abstractproperty
    def pseudo(self):

    @abc.abstractproperty
    def generator_code(self):

    @abc.abstractmethod
    def get_code_imput(self):
        """

        """






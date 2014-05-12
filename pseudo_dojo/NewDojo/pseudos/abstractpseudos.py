"""
here we define the main abstract class for a pp and the main derived abstact pp's
the concrete pp's should be stored with their data in subfolderes
"""

__author__ = 'setten'

from abc import ABCMeta, abstractproperty, abstractmethod
from pseudo_dojo.NewDojo.codes.generatorcodes import get_generator_interface
from pseudo_dojo.NewDojo.codes.systemcodes import get_system_code_interface


class AbstractPP(object):
    """
    Main abstract pp class defining the interface
    only the methods defined here are to be used in the dojo
    """
    __metaclass__ = ABCMeta

    def __init__(self):
        self.generator_interface = get_generator_interface(self.generator_code)
        self._data_set = None

    def create_data_set(self):
        """
        method to create by code the data_set from definition
        """
        self._data_set = self.generator_interface.create_data_set(self.definition)

    def write_pp_file(self, path, system_code):
        """
        method to write the pp to the give path to be used by the system code
        """
        code_interface = get_system_code_interface(system_code)
        code_interface.write_pseudo(self.data_set, path)

    @property
    def data_set(self):
        return self._data_set

    @abstractproperty
    def definition(self):
        """
        the set of parameters that define the pp, passing these to the defined code should suffice to generate the pp
        """

    @abstractproperty
    def generator_code(self):
        """
        property that contains the name of the generator code
        """


class AbstractDefinition(object):
    """

    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def read_from_file(self, path):
        """
        method to read the defining parameters from a file
        """

    @abstractmethod
    def write_fo_file(self, path):
        """
        method to write the defining parameters to file
        """


class NCDefinition(AbstractDefinition):
    """

    """
    def __init__(self):
        self.generator_code = 'some code'
        """
        here we should define all the properties defining a NC
        """

    def read_from_file(self, path):
        raise NotImplementedError

    def write_fo_file(self, path):
        raise NotImplementedError


class AbstractDataSet(object):
    """

    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def read_from_file(self, path):
        """
        method to read an actual data_set from file
        """

    @abstractmethod
    def write_to_file(self, path):
        """
        method to write an actual data_set from file
        """


class NCDataSet(AbstractDataSet):
    """

    """

    def read_from_file(self, path):
        """
        here we need to implement how to read an NCDataSet from file
        """
        raise NotImplementedError

    def write_to_file(self, path):
        """
        here we need to implement how to write a NCDataSet to file
        """
        raise NotImplementedError


class USDataSet(AbstractDataSet):
    """

    """

    def read_from_file(self, path):
        raise NotImplementedError

    def write_to_file(self, path):
        raise NotImplementedError


class PAWDataSet(AbstractDataSet):
    """

    """

    def read_from_file(self, path):
        raise NotImplementedError

    def write_to_file(self, path):
        raise NotImplementedError


class GeneralNCPP(AbstractPP):
    """
    prototype of specific pp of NC type
    used by the factory functions and as a template for explicitly programmed NCPP's
    """
    def __init__(self):
        super(self.__class__, self).__init__()
        self._data_set = NCDataSet()
        self._definition = NCDefinition()                            # the definition may become public at some point
        self._code = 'name'

    @property
    def definition(self):
        return self._definition

    @property
    def generator_code(self):
        return self._code


######
# API
######
# factory functions


def get_ncpp_from_definition_file(path):
    """
    factory function to create an NCPP from a file containing the definition
    """
    ncpp = GeneralNCPP()
    ncpp.definition.read_from_file(path)
    ncpp._code = ncpp.definition.generator_code
    ncpp.create_data_set()
    return ncpp


def get_ncpp_from_data_set_file(path):
    """
    factory function to create an NCPP from a file containing a data_set
    the definition is not know in this wa
    """
    ncpp = GeneralNCPP()
    ncpp.data_set.read_from_file(path)
    ncpp._code = ncpp.data_set.generator_code
    ncpp._definition = None
    return ncpp

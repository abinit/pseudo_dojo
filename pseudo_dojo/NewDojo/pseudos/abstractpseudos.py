"""
here we define the main abstract class for a pp and the main derived abstact pp's
the concrete pp's should be stored with their data in subfolderes
"""

__author__ = 'setten'

from abc import ABCMeta, abstractproperty, abstractmethod
from pseudo_dojo.NewDojo.codes.generatorcodes import get_generator_interface


class AbstractPP(object):
    """
    Main abstract pp class defining the interface
    only the methods defined here are to be used in the dojo
    """
    __metaclass__ = ABCMeta

    @abstractproperty
    def definition(self):
        """
        the set of parameters that define the pp, passing these to the defined code should suffice to generate the pp
        """

    @abstractproperty
    def data_set(self):
        """
        the actual pp
        """

    @abstractproperty
    def code(self):
        """
        property that contains the name of the generator code
        """

    @abstractmethod
    def read(self):
        """
        method to read a pp from a database and put it in self.dataset
        """

    @abstractmethod
    def write(self, path):
        """
        method to write the pp to the give path to be used by the system code
        """

    @abstractmethod
    def create(self):
        """
        method to create by code the data_set from definition
        """


class AbstractDefinition(object):
    """

    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def read_definition_from_file(self, path):
        """
        method to read the defining paramers from a file
        """


class NCDefinition(AbstractDefinition):
    """

    """
    def __init__(self):
        """

        """

    def read_definition_from_file(self, path):
        pass


class AbstractDataSet(object):
    """

    """
    __metaclass__ = ABCMeta


class NCDataSet(AbstractDataSet):
    """

    """


class ANCPP(AbstractPP):
    """
    prototype of specific pp of NC type
    """
    def __init__(self):
        self._data_set = NCDataSet()
        self._definition = NCDefinition()                            # the definition may become public at some point
        self._code = 'name'
        self.creator_interface = get_generator_interface(self.code)
        # now we should define the pp here by
        #   reading a definition
        #   or hardcoding a definition
        #   and creating the dataset
        # or
        #   read a data_set from file

    @property
    def definition(self):
        return self._definition

    @property
    def data_set(self):
        return self._data_set

    @property
    def code(self):
        return self._code

    def create(self):
        self._data_set = self.creator_interface.create(self.definition)

    def write(self, path):
        pass

    def read(self):
        pass
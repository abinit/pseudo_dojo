__author__ = 'setten'

from pseudo_dojo.NewDojo.pseudos.abstractpseudos import NCDataSet, NCDefinition, AbstractPP
from pseudo_dojo.NewDojo.codes.generatorcodes import get_generator_interface


class SomeNCPP(AbstractPP):
    """
    prototype of specific pp of NC type
    used as a template for explicitly programmed NCPP's
    """
    def __init__(self):
        super(self.__class__, self).__init__()
        self._data_set = NCDataSet()
        self._definition = NCDefinition()                            # the definition may become public at some point
        self._code = 'name'
        # now we should define the pp here by
        #   reading a definition:
        #       self.definition.read_from_file('path')
        #   or hardcoding a definition
        #   and creating the data_set:
        #       self.create()
        # or
        #   read a data_set from file
        #       self._data_set.read_from_file('path')

    @property
    def definition(self):
        return self._definition

    @property
    def generator_code(self):
        return self._code
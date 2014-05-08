__author__ = 'setten'

from pseudo_dojo.NewDojo.pseudos.abstractpseudos import NCDataSet, NCDefinition, AbstractPP
from pseudo_dojo.NewDojo.codes.generatorcodes import get_generator_interface


class SomeCPP(AbstractPP):
    """
    prototype of specific pp of NC type
    the doc string should contain a human readable concise description
    """
    def __init__(self):
        self._data_set = NCDataSet()
        self._definition = NCDefinition()
        self._code = 'name'
        self.creator_interface = get_generator_interface(self.generator_code)

    @property
    def definition(self):
        return self._definition

    @property
    def data_set(self):
        return self._data_set

    @property
    def generator_code(self):
        return self._code

    def create(self):
        self._data_set = self.creator_interface.create(self.definition)

    def write(self, path):
        pass

    def read(self):
        pass
__author__ = 'setten'

from abc import ABCMeta, abstractproperty, abstractmethod
from pseudo_dojo.NewDojo.codes.systemcodes import get_system_code_interface, AbstractSystemCodeInterface


class AbstractKata(object):
    """
    Abstract class for defining the interface for kata's (tests)
    the dojo should only use this interface in implementing
    """
    __Metaclass__ = ABCMeta

    def set_code_interface(self, code):
        self.code_interface = get_system_code_interface(code)

    def test_code_type(self):
        if isinstance(super(self.code_interface), self.code_interface):
            return True

    @abstractproperty
    def code_type(self):
        """
        the type of code that is needed by this code given by the abstact base class of that type i.e.
        AbstractSystemCode() of AbstractAtomicCode()
        """

    @abstractmethod
    def hajime(self):
        """
        run the test
        """


class DeltaFactorKata(AbstractKata):
    """

    """

    @property
    def code_type(self):
        return AbstractSystemCodeInterface

    def hajime(self):
        """

        """


######
# API
######


def get_kata(kata):
    """
    factory function to return an instance of a systems code interface
    """
    cls = {'deltafactorkata': DeltaFactorKata()}
    return cls[kata]()
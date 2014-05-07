__author__ = 'setten'

from pseudo_dojo.NewDojo.katas.abstractkata import get_kata
from pseudo_dojo.NewDojo.codes.systemcodes import get_system_code_interface


class Dojo(object):
    """
    class for combining a set of shiai's for differend kata's, codes and pp's
    """


class BaseShiai(object):
    """
    (match, competition), object that pacts all together for the excecution of a specific test:

    kata : the actual definition of the test
    code : the code that is to be used
    pp   : the pp that is to be tested

    this is the base class for shiai's creating the general framework
    from this one we can extend

    """
    def __init__(self, kata, code, pp):
        # these are just string names
        self._code = code
        self._kata = kata
        self._pp = pp

        # make self.kata an instance of the kata we want
        self.kata = get_kata(self._kata)
        # give the kata the code interface
        self.kata.set_code_interface(self._code)

        # todo pass the pp in some way to the code
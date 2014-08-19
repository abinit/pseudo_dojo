__author__ = 'setten'

from time import sleep
from pseudo_dojo.NewDojo.katas.abstractkata import get_kata
from pseudo_dojo.NewDojo.codes.systemcodes import get_system_code_interface


class DataBase(object):
    """
    class for an active (FW) data_base
    """

    def initialize(self):
        """
        initialize a FW database that is going to be treated as a subject
        """

    def add_work(self, work):
        """
        add a work to the database
        """

    def has_new_data(self, dojo_id):
        """
        check if the data_base contains new data since the last call
        """

    def get_new_dojo_id(self):
        """
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
        self.data = {}

        # make self.kata an instance of the kata we want
        # for the moment have only single element kata's
        self.kata = get_kata(self._kata, pp)
        # give the kata the code interface
        self.kata.set_code_interface(self._code)

        # todo pass the pp in some way to the code


class Dojo(object):
    """
    class for combining a set of shiai's for differend kata's, codes and pp's
    """

    def __init__(self):
        self.data_base = DataBase()
        self.katas = []
        self.codes = ['abinit']
        self.pps = []
        self.shiais = {}

    def update_katas(self):
        """
        method to updata the kata's to be performed
        """

    def update_codes(self):
        """
        method to update the codes to be used
        """

    def update(self):
        self.update_codes()
        self.update_katas()

    def simple_collection_run(self):
        self.data_base.initialize()
        self.update()
        all_done = False

        # submit all calculations to the database:
        for code in self.codes:
            for pp in self.pps:
                for kata in self.katas:
                    # for the moment we only have only single element katas
                    shiai = BaseShiai(kata=kata, code=code, pp=pp)
                    self.data_base.add_work(shiai, s_id)
                    self.shiais.update({'s_id': s_id, 'shiai': shiai, 'done': False})

        # start watching the data_base for complete data sets
        while not all_done:
            if self.data_base.has_new_data():  # todo find out how to do this ...
                for shiai in self.shiais:
                    pass
                    # when it is complete and not done
                    # evaluate the data
            else:
                sleep(5)


if __name__ == "__main__":
    my_dojo = Dojo()
    my_dojo.simple_collection_run()

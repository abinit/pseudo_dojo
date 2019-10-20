#!/usr/bin/env python
import os.path
import sys

from pprint import pprint

import nist_database


def sort_dict(d, key=None, reverse=False):
    """
    Returns an OrderedDict object whose keys are ordered according to their
    value.

    Args:
        d:
            Input dictionary
        key:
            function which takes an tuple (key, object) and returns a value to
            compare and sort by. By default, the function compares the values
             of the dict i.e. key = lambda t : t[1]
        reverse:
            allows to reverse sort order.
    """
    import collections
    kv_items = [kv for kv in d.items()]

    # Sort kv_items according to key.
    if key is None:
        kv_items.sort(key=lambda t: t[1], reverse=reverse)
    else:
        kv_items.sort(key=key, reverse=reverse)

    # Build ordered dict.
    return collections.OrderedDict(kv_items)

#class AtomicData(dict):
#
#     _mandatory_keys = [
#         "Etot",
#         "Ekin",
#         "Ecoul",
#         "Eenuc",
#         "Exc",
#     ]
#
#     def __init__(self, *args, **kwargs):
#        super(AtomicData, self).__init__(*args, **kwargs)
#
#        for k in self._mandatory_keys:
#            if k not in self:
#                raise ValueError("mandatory key %s is missing" % k)
#        #self.symbol
#        #self.Z = Z
#        #self.configurations = configurations


def make_nist_configurations(path):

    neutral, cations, symb2Z = parse_nist_configurations(path)

    print("# Computer generated code")
    print("neutral = ")
    pprint(neutral)
    print("")

    print("# Computer generated code")
    print("cations = ")
    pprint(cations)
    print("")

    print("# Computer generated code")
    print("symb2Z = ")
    pprint(symb2Z)
    print("")


def parse_nist_configurations(path):
    """Read and parse the file with the configurations."""
    # Z    Symbol     Neutral                    Positive ion
    # 1      H       1s^1                           -
    # 2      He      1s^2                        1s^1
    neutral, cations, symb2Z = {}, {}, {}

    count = 1
    with open(os.path.join(path, "configurations"), "r") as fh:
        for line in fh:
            if not (len(line) > 1 and line[1].isdigit()):
                continue
            head, catio = line[:44], line[44:].strip()
            toks = head.split()
            Z, symbol, neutr = int(toks[0]), toks[1], " ".join(t for t in toks[2:])
            assert Z == count
            count += 1

            neutr = neutr.replace("^", "")
            catio = catio.replace("^", "")

            neutral[symbol] = neutr
            cations[symbol] = catio if catio != "-" else None
            symb2Z[symbol] = Z

        return neutral, cations, symb2Z


def occupations_from_symbol(symbol):
    Z = symb2Z[symbol]
    configuration = neutral[Z]

    if configuration[0][0] == '[':
        occupations = occupations_from_symbol(configuration[0][1:-1])
        configuration = configuration[1:]
    else:
        occupations = []
    for s in configuration:
        occupations.append((s[:2], int(s[3:])))
    return occupations


def extract_nistdata(path, nist_xctype, iontype):
    # http://www.physics.nist.gov/PhysRefData/DFTdata/

    # Read and parse the configurations file
    # Z    Symbol     Neutral                    Positive ion
    # 1      H       1s^1                           -
    # 2      He      1s^2                        1s^1
    Ztable = {}
    configurations = [['X', '']]
    Z = 1

    with open(os.path.join(path, "configurations"), "r") as fh:
        for line in fh:
            if len(line) > 1 and line[1].isdigit():
                line = line[:44].split()
                symbol = line[1]
                Ztable[symbol] = Z
                assert int(line[0]) == Z
                configurations.append(line[2:])
                Z += 1

    def occupations_from_symbol(symbol):
        Z = Ztable[symbol]
        configuration = configurations[Z]
        if configuration[0][0] == '[':
            occupations = occupations_from_symbol(configuration[0][1:-1])
            configuration = configuration[1:]
        else:
            occupations = []
        for s in configuration:
            occupations.append((s[:2], int(s[3:])))
        return occupations

    # Ex of file with AE results:
    # Etot  =     -675.742283
    # Ekin  =      674.657334
    # Ecoul =      285.206130
    # Eenuc =    -1601.463209
    # Exc   =      -34.142538
    # 1s       -143.935181
    # 2s        -15.046905
    # 2p        -12.285376
    # 3s         -1.706331
    # 3p         -1.030572
    # 4s         -0.141411
    nistdata = {}
    spdf = {'s': 0, 'p': 1, 'd': 2, 'f': 3}

    for (symbol, Z) in Ztable.items():
        #print("in symbol %s" % symbol)
        occupations = occupations_from_symbol(symbol)

        fname = os.path.join(path, nist_xctype, iontype, '%02d%s' % (Z, symbol))

        energies, eigens = {}, {}

        with open(fname, 'r') as fh:
            for n in range(5):
                ename, evalue = fh.readline().split("=")
                energies[ename.strip()] = float(evalue)

            for line in fh:
                state, eig = line.split()
                eigens[state] = float(eig)

        nloe = []
        for (state, occ) in occupations:
            n = int(state[0])
            l = spdf[state[1]]
            eig = eigens[state]
            nloe.append((n, l, occ, eig))

        nistdata[symbol] = (Z, nloe, energies)

    return nistdata


def make_nistmodule(path, nist_xctype):
    #for iontype in ["neutrals", "cations",]:
    for iontype in ["neutrals",]:
        print('# Computer generated code: nist_xctype = %s, iontype = %s' % (nist_xctype, iontype))
        print("format:\n\t element: (atomic number, [(n, l, occ, energy), ...]")
        print('%s = ' % iontype)
        data = extract_nistdata(path, nist_xctype, iontype)
        pprint(data)
        #print(json.dumps(data, indent=4))


if __name__ == '__main__':
    import sys
    from qatom import states_from_string, AtomicConfiguration

    for symbol in nist_database.symbols:
        aconf = AtomicConfiguration.neutral_from_symbol(symbol)
        print(aconf)
    sys.exit(1)

    path = sys.argv[1]

    #neutral = {}
    #for (symbol, confstr) in nist_database.neutral.items():
    #    states = states_from_string(confstr)
    #    print(confstr)
    #    neutral[symbol] = (nist_database.symb2Z[symbol], states)
        #for s in states:
        #    pprint(tuple(s))
        #aconf = AtomicConfiguration.from_string(confstr, symb2Z[symbol])

    #make_nist_configurations(path)

    #xctypes = ["LDA", "ScRLDA"] #, "LSD", "RLDA",]
    nist_xctypes = ["LDA",]
    for xctype in nist_xctypes:
        print("xctype: %s" % xctype)
        make_nistmodule(path, xctype)

    #iontype = "neutrals"
    #dft_data = extract_nistdata(sys.argv[1], nist_xctype, iontype)
    #make_nistmodule(sys.argv[1], "LDA")

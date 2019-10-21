#!/usr/bin/env python


class AtomicData:
    def __init__(self, name, Z, mass, radius, configuration):
        self.name = name
        self.Z = Z
        self.mass = mass
        self.radius = radius
        self.configuration = configuration


def core_states(symbol):
    """Method returning the number of core states for given element."""
    core = parameters[symbol].get('core', '')

    # Parse core string:
    j = 0
    if core.startswith('['):
        a, core = core.split(']')
        core_symbol = a[1:]
        j = len(configurations[core_symbol][1])

    Njcore = j + len(core) // 2
    return Njcore


if __name__ == '__main__':
    import pprint
    import sys
    import os.path
    # http://www.physics.nist.gov/PhysRefData/DFTdata/
    path = os.path.abspath(sys.argv[1])

    Ztable = {}
    configurations = [['X', '']]
    Z = 1
    for line in open(os.path.join(path, 'configurations')):
        if len(line) > 1 and line[1].isdigit():
            line = line[:44].split()
            symbol = line[1]
            Ztable[symbol] = Z
            configurations.append(line[2:])
            Z += 1

    def get_occupations(symbol):
        Z = Ztable[symbol]
        configuration = configurations[Z]
        if configuration[0][0] == '[':
            occupations = get_occupations(configuration[0][1:-1])
            configuration = configuration[1:]
        else:
            occupations = []
        for s in configuration:
            occupations.append((s[:2], int(s[3:])))
        return occupations

    dftdata = {}
    spdf = {'s': 0, 'p': 1, 'd': 2, 'f': 3}
    for symbol, Z in Ztable.items():
        occupations = get_occupations(symbol)
        f = open(os.path.join(path, 'LDA', 'neutrals', '%02d%s' % (Z, symbol)), 'r')
        for n in range(5):
            f.readline()
        epsilons = {}
        for line in f:
            state, eps = line.split()
            epsilons[state] = float(eps)
        nloe = []
        for state, occ in occupations:
            n = int(state[0])
            l = spdf[state[1]]
            eps = epsilons[state]
            nloe.append((n, l, occ, eps))
        dftdata[symbol] = (Z, nloe)
    print('# Computer generated code:')
    print("")
    print('configurations = ')
    pprint.pprint(dftdata)
    print("")

from __future__ import division, print_function

import os
import sys
import time
import collections
import shutil
import json
import numpy as np

from pprint import pprint
from subprocess import Popen, PIPE

from pseudo_dojo.core import QState, AtomicConfiguration, plot_logders
from pseudo_dojo.ppcodes.ape.apeio import (ape_read_waves, ape_read_potentials, ape_read_densities, 
                                           ape_read_logders, ape_read_dipoles, ape_check_ppeigen, ape_check_ghosts)

__version__ = "0.1"

##########################################################################################

class ApeVariable(object):

    def __init__(self, name, type, doc, allowed_values=None):
        self.name = name
        self.type = type
        self.doc = doc
        self.allowed_values = allowed_values if allowed_values is not None else []

    def set_value(self, value):
        self._value = value

    @property
    def value(self):
        try:
            return self._value
        except AttributeError:
            return None

    def isvalid(self):
        if self.allowed_values:
            return self.value in self.allowed_values
        else:
            return True

    def to_apeinput(self):
        return "%s = %s" % (self.name, self.value)


class ApeBlock(object):

    def __init__(self, name, types, doc, allowed_values=None):
        self.name = name
        self.types = types
        self.doc = doc
        self.allowed_values = allowed_values if allowed_values is not None else []

        self._table = []

    def add_entry(self, entry):
        row = []
        for i, item in enumerate(entry):
            row.append(self.types[i](item))
        self._table.append(row)

    def to_apeinput(self):
        lines = ["\%%s" % self.name]
        for (state, rcut, scheme) in zip(self.states, self.core_radii, self.schemes):
            lines += ["%s | %s | %s | %s" % (state.n, state.l, rcut, scheme)]
        lines += ["%"]
        return "\n".join(lines)

##########################################################################################

class ApeAtomicConfiguration(AtomicConfiguration):

    def to_apeinput(self):
        lines  = ["NuclearCharge = %s " % self.Z]
        lines += ["SpinMode = %s" % self.spin_mode]
        lines += ["%Orbitals"]
        for state in self:
            lines += state.to_apeinput()
        lines += ["%"]
        return lines

    @classmethod
    def from_apeinput(cls, lines):

        if isinstance(lines, str):
            with open(lines, "r") as fh:
                lines = fh.readlines()
        # TODO
        # Example (spin unpolarized, no spinor)
        #%Orbitals
        #1 | 0 | 2.0 
        #2 | 0 | 2.0 
        #2 | 1 | 6.0 
        #3 | 0 | 2.0 
        #3 | 1 | 2.0 
        #%
        START_ORBS = "%Orbitals"
        for (lineno, line) in enumerate(lines):
            if line.startswith("NuclearCharge"): s, Z = line.split("=")
            if line.startswith(START_ORBS): start = lineno

        ape_states, in_orbitals = [], False
        for line in lines[start+1:]:
            line = line.strip()
            if line.startswith("%"): break
            ape_states.append(line)
        assert ape_states

        states = []
        for (i, s) in enumerate(ape_states):
            print("s:",s)

            if "|" not in s:
                # Handle noble gas.
                assert i == 0
                s = s.translate(None, "'\"")
                noble_gas = ApeAtomicConfiguration.neutral_from_symbol(s)
                states = noble_gas.states
            else:
                toks = s.split("|")
                assert len(toks) == 3
                qstate = QState(n=toks[0], l=toks[1], occ=toks[2])
                states.append(qstate)

        return cls(Z, states)

##########################################################################################

class ApeRadialMesh(dict):
    """
    The radial mesh used by APE to represent radial functions. 
    """
    # Supported variables
    _VARS = [
        "MeshType",
        "MeshStartingPoint",
        "MeshOutmostPoint",
        "MeshNumberOfPoints", 
        "MeshDerivMethod",
        "MeshFiniteDiffOrder",
    ]

    def __init__(self, *args, **kwargs):
        super(ApeRadialMesh, self).__init__(*args, **kwargs)

        for k in self:
            if k not in self._VARS:
                raise ValueError("%s is not a registered key" % k)

    def to_apeinput(self):
        return ["%s = %s" % kv for kv in self.items()]


class ApeControl(dict):
    """
    The self consistent field procedure will stop when one of the convergence criteria is fulfilled. 
    At each iteration the new guess potential is built mixing the input and output potentials.
    """
    # Supported variables
    _VARS = [
        # SCF
        "MaximumIter", 
        "ConvAbsDens",
        "ConvRelDens",
        "ConvAbsEnergy",
        "ConvRelEnergy",
        "SmearingFunction",   # This one seems a very delicate point
        "MixingScheme",
        "Mixing",
        "MixNumberSteps",
        "MaximumIter",
        # Eigensolver
        "EigenSolverTolerance",
        "ODEIntTolerance",
        "ODESteppingFunction",
        "ODEMaxSteps",
        "EigenSolverTolerance",
        # Generic
        "Verbose",
    ]

    def __init__(self, **kwargs):
        super(ApeControl, self).__init__(**kwargs)
                                                                   
        for k in self:
            if k not in self._VARS:
                raise ValueError("%s is not a registered key" % k)
                                                                   
    def to_apeinput(self):
        return ["%s = %s" % kv for kv in self.items()]



class ApePPComponents(object):

    @classmethod
    def from_strings(cls, sep="|", *strings):
        """
        Instanciate the object from a list of strings 
        Example: "3s:1.2:tm"
        """
        states, core_radii, schemes = [], [], []

        occ0 = 0.0
        for s in strings:
            try:
                tokens = s.split(sep)
                assert len(tokens) == 3
                assert len(tokens[0]) == 2
                n = tokens[0][0]
                l = tokens[0][1]
                rc = float(tokens[1])
                scheme = tokens[2]

                # Add them to the lists.
                states.append(QState(n, l, occ0))
                core_radii.append(rc)
                schemes.append(scheme)

            except:
                raise ValueError("Malformatted string %s" % s)

        return cls(states, core_radii, schemes)

    def __init__(self, states, core_radii, schemes):
        self.states = states
        self.core_radii = core_radii
        self.schemes = schemes

    def to_apeinput(self):
        lines = ["%PPComponents"]
        for (state, rcut, scheme) in zip(self.states, self.core_radii, self.schemes):
            lines += ["%s | %s | %s | %s" % (state.n, state.l, rcut, scheme)]
        lines += ["%"]
        return lines

##########################################################################################

class ApePPSetup(object):

    def __init__(self, pp_components, core_correction=0, llocal=-1):
        self.pp_components = pp_components
        self.core_correction = core_correction 
        self.llocal = llocal

    def to_apeinput(self):
        lines =  ["# PseudoPotential Setup"]
        lines += ["CoreCorrection = %s" % self.core_correction]
        lines += ["Llocal = %s" % self.llocal]
        lines += self.pp_components.to_input()
        return lines

##########################################################################################

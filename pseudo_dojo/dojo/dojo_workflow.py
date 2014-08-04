"""Base class for Dojo Workkflows."""
from __future__ import division, print_function

import abc

from pymatgen.io.abinitio.workflows import Workflow

class DojoWorkflow(Workflow):
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def pseudo(self):
        """Pseudo"""

    @abc.abstractproperty
    def dojo_trial(self):
        """String identifying the DOJO trial. Used to write results in the DOJO_REPORT."""

    @property
    def dojo_accuracy(self):
        return self._dojo_accuracy

    def set_dojo_accuracy(self, accuracy):
        self._dojo_accuracy = accuracy

    def write_dojo_report(self, report, overwrite_data=False, ignore_errors=False):
        """Write/update the DOJO_REPORT section of the pseudopotential."""
        # Read old_report from pseudo.
        old_report = self.pseudo.read_dojo_report()

        dojo_trial, dojo_accuracy = self.dojo_trial, self.dojo_accuracy
        if dojo_trial not in old_report:
            # Create new entry
            old_report[dojo_trial] = {}
        else:
            # Check that we are not going to overwrite data.
            if self.dojo_accuracy in old_report[dojo_trial] and not overwrite_data:
                raise RuntimeError("%s already exists in DOJO_REPORT. Cannot overwrite data" % dojo_trial)

        # Update old report card with the new one.
        old_report[dojo_trial][dojo_accuracy] = report

        # Write new report
        self.pseudo.write_dojo_report(old_report)

# coding: utf-8
"""Test PseudoDojo command line scripts."""

import os
import pseudo_dojo.data as pdj_data

from scripttest import TestFileEnvironment
from monty.inspect import all_subclasses
from pseudo_dojo.core.testing import PseudoDojoTest
from pseudo_dojo.pseudos import dojotable_absdir

script_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))


def test_if_all_scripts_are_tested():
    """Testing if all scripts are tested"""
    tested_scripts = set(os.path.basename(c.script) for c in all_subclasses(ScriptTest))
    all_scripts = set(f for f in os.listdir(script_dir) if f.endswith(".py"))
    not_tested = all_scripts.difference(tested_scripts)

    if not_tested:
        print("The following scripts are not tested")
        for i, s in enumerate(not_tested):
            print("[%d] %s" % (i, s))

    #assert not_tested == set(['dojogbrv.py'])


class ScriptTest(PseudoDojoTest):
    loglevel = "--loglevel=ERROR"
    verbose = "--verbose"

    def get_env(self):
        #import tempfile
        #env = TestFileEnvironment(tempfile.mkdtemp(suffix='', prefix='test_' + script))
        env = TestFileEnvironment()

        # Use Agg backend for plots.
        #env.writefile("matplotlibrc", "backend : Agg")

        # Start with --help. If this does not work...
        env.run(self.script, "--help")

        # Script must provide a version option
        #r = env.run(self.script, "--version", expect_stderr=True)
        #assert r.stderr.strip() == "%s version %s" % (os.path.basename(self.script), abilab.__version__)
        return env


class TestDojoRun(ScriptTest):
    script = os.path.join(script_dir, "dojorun.py")

    def test_dojorun_dryrun(self):
        """Testing dojorun.py script in dry-run mode."""
        env = self.get_env()

        # Build new env to run the script in dry-run mode
        # Copy file from data to env.
        env = TestFileEnvironment(template_path=pdj_data.dirpath)
        env.writefile("Si.psp8", frompath="Si.psp8")
        env.writefile("Si.djrepo", frompath="Si.djrepo_empty")
        env.run(self.script, "Si.psp8", self.loglevel, self.verbose, "--dry-run")


class TestDojoData(ScriptTest):
    script = os.path.join(script_dir, "dojodata.py")

    def test_dojodata(self):
        """Testing dojodata.py script"""
        env = self.get_env()

        # Test dojodata table
        path = dojotable_absdir("ONCVPSP-PBE-PDv0.3")
        env.run(self.script, "table", path, self.loglevel, self.verbose)

        # Test dojodata check
        #path = dojotable_absdir("ONCVPSP-PBE-PDv0.3")
        #env.run(self.script, "check", path, self.loglevel, self.verbose)


class TestOncv(ScriptTest):
    script = os.path.join(script_dir, "dojoncv.py")

    def test_oncv(self):
        """Testing dojodata.py script"""
        env = self.get_env()
        #env.run(self.script, self.loglevel, "man", "ecut")


class TestDojoGbrv(ScriptTest):
    script = os.path.join(script_dir, "dojogbrv.py")

    def test_dojogbrv_dryrun(self):
        """Testing dojogbrv.py script."""
        env = self.get_env()

        env.run(self.script, "info", self.loglevel, self.verbose)
        env.run(self.script, "find", "Si", "O", self.loglevel, self.verbose)

        # TODO
        # Build new env to run the script in dry-run mode
        # Copy file from data to env.
        #env = TestFileEnvironment(template_path=pdj_data.dirpath)
        #env.writefile("Si.psp8", frompath="Si.psp8")
        #env.writefile("Si.djrepo", frompath="Si.djrepo_empty")
        #env.run(self.script, "runps", "Si.psp8", self.loglevel, self.verbose, "--dry-run")

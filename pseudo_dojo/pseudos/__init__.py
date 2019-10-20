"""
Functions providing access to file data. Mainly used to build the public APIs and write unit tests.
"""
import os

from monty.termcolor import cprint
from pseudo_dojo.core.pseudos import dojopseudo_from_file

DOJOTABLE_BASEDIRS = [
    "ONCVPSP-PBE-PDv0.2",
    "ONCVPSP-PBE-PDv0.3",
    "ONCVPSP-PBEsol-PDv0.3",
    "ONCVPSP-PW-DEV",
    # Tables for the paper
    "ONCVPSP-PBE-PDv0.4",
    "ONCVPSP-PW-PDv0.4",
    "ONCVPSP-PBEsol-PDv0.4"
]

here = os.path.dirname(__file__)


def dojotable_absdir(basedir):
    """
    Return the absolute dirpath of the table from its basename
    """
    if basedir not in DOJOTABLE_BASEDIRS:
        raise RuntimeError(
           "%s is not registered in DOJOTABLE_BASEDIRS\n"
           "Change pseudo_dojo/pseudos/__init__.py" % basedir)

    return os.path.join(here, basedir)


def all_dojotable_absdirs():
    """
    List with the absolute path of the directories containing `stable` pseudos"
    """
    return [dojotable_absdir(bdir) for bdir in DOJOTABLE_BASEDIRS]


def as_dojo_path(path):
    """
    Translate an absolute path into an official path inside the pseudo_dojo package.

    This function is extremaly useful when we have absolute paths
    produced on a different machine and we want to get the corresponding
    path on the local host.
    """
    if os.path.exists(path): return path
    # Use the final part of the path and change the root.
    head, base = os.path.split(path)
    _, dirname = os.path.split(head)
    dojodir = dojotable_absdir(dirname)
    path = os.path.join(dojodir, base)
    if not os.path.exists(path):
        raise RuntimeError("After filepath mangling. No such file: %s" % path)
    return path


def check_pseudo_path(path, verbose=0):
    """
    Check a pseudopotential given the filepath. Warnings are printed to stdout.
    Return 0 if success.
    """
    pseudo = dojopseudo_from_file(path)
    if pseudo is None:
        cprint("[%s] Pseudo.from_file returned None. Something wrong in file!" % path, "red")
        return 1

    return check_pseudo(pseudo, verbose=verbose)


def check_pseudo(pseudo, check_trials=None, verbose=0):
    """
    Check a pseudopotential object. Warnings are printed to stdout
    Return 0 if success.
    """
    retcode = 0
    try:
        report = pseudo.dojo_report
    except Exception as exc:
        cprint("Connot find dojo_report associated to: [%s]" % os.path.relpath(pseudo.filepath), "red")
        if verbose: print("Python Exception:\n%s", exc)
        retcode += 1
        return retcode

    if "ppgen_hints" not in report:
        cprint("[%s] old version without ppgen_hints" % os.path.relpath(pseudo.filepath), "red")
        retcode += 1

    if report["version"] != "1.0":
        cprint("[%s] wrong version: %s" % (os.path.relpath(pseudo.filepath), report["version"]), "red")
        retcode += 1

    if report["md5"] != pseudo.compute_md5():
        cprint("Incosistent md5 in [%s]" % os.path.relpath(pseudo.filepath), "red")
        retcode += 1

    if report["symbol"] != pseudo.symbol:
        cprint("Inconsistent symbol in [%s]" % os.path.relpath(pseudo.filepath), "red")
        retcode += 1

    # This part is commented because we are gonna refactor the DojoReport
    try:
        error = report.check(check_trials=check_trials)
        #error = report.check(check_trials=["deltafactor"])
        #error = report.check(check_trials=["deltafactor", "gbrv_bcc", "gbrv_fcc"])
        if error:
            retcode += 1
            cprint("Invalid DojoReport in [%s]" % os.path.relpath(pseudo.filepath), "red")
            if verbose: print(error)

    except Exception as exc:
        retcode += 1
        cprint("Python exception in [%s]" % os.path.relpath(pseudo.filepath), "red")
        if verbose: print(str(exc))

    if retcode != 0:
        cprint("[%s] is not valid" % os.path.relpath(pseudo.filepath), "red")
        if not verbose: cprint("Use --verbose for more info", "yellow")

    return retcode


#def check_djrepo(djrepo, verbose=0):
#    """
#    Check a pseudopotential given the filepath of the djrepo file.
#    Warnings are printed to stdout.
#    Return 0 if success.
#    """
#    # Get the path of the pseudo from djrepo
#    import json
#    djrepo = os.path.abspath(djrepo)
#    with open(djrepo, "rt") as fh:
#        d = json.load(fh)
#        path = os.path.join(os.path.dirname(djrepo), d["basename"])
#
#    return check_pseudo_path(path, verbose=verbose)

#!/usr/bin/env python
"""Script used by maintainers to check djrepo files."""
import sys
import os
import argparse
import glob

from monty.functools import prof_main
from monty.termcolor import cprint
from monty.os.path import find_exts
from pseudo_dojo.core.pseudos import dojopseudo_from_file
from pseudo_dojo.core.dojoreport import DojoReport, dojo_dfact_results


def djrepo_fix(options):
    """
    This is a maintentance tool to:

        #. Regenerate the md5 value in the DojoReport file.
        #. Replace norm-conserving with NC.
    """
    def fix_filepath(pp_filepath):
        pseudo = dojopseudo_from_file(pp_filepath)
        if pseudo is None:
            print("Error while parsing %s" % pp_filepath)
            return 1

        # Change md5
        pseudo.dojo_report["md5"] = pseudo.compute_md5()
        if pseudo.dojo_report["pseudo_type"] == "norm-conserving":
            pseudo.dojo_report["pseudo_type"] = "NC"

        # Add basename
        if "basename" not in pseudo.dojo_report:
            pseudo.dojo_report["basename"] = pseudo.basename
        if pseudo.dojo_report["basename"] != pseudo.basename:
            print("Inconsistent basename in %s" % pp_filepath)
            return 1

        # Remove ebands
        #pseudo.dojo_report.pop("ebands", None)
        #if pseudo.dojo_report.pop("phonon", None) is None: return 0

        # Rename entries in FR pseudos.
        #oldnew = [("deltafactor", "deltafactor_soc"),
        #          ("gbrv_fcc", "gbrv_fcc_soc"),
        #          ("gbrv_bcc", "gbrv_bcc_soc")]
        #for old, new in oldnew:
        #    if old in pseudo.dojo_report:
        #        pseudo.dojo_report[new] = pseudo.dojo_report.pop(old)
        #pseudo.dojo_report.reorder()

        pseudo.dojo_report.json_write()
        return 0

    print("Will fix %s djrepo files" % len(options.paths))
    retcode = 0
    for path in options.paths:
        if options.verbose: print(path)
        retcode += fix_filepath(path)
        if retcode:
            cprint("retcode != 0 in %s" % path, "red")
            break

    return retcode


def djrepo_check(options):
    """Check djrepo files."""
    retcode = 0
    for path in options.paths:
        #print(path)
        pseudo = dojopseudo_from_file(path)
        report = pseudo.dojo_report
        validation = report.get("validation", None)
        if validation is None:
            retcode += 1
            cprint("[%s]: no validation entry" % path)
            continue
        hints = report.get("hints", None)
        if hints is None:
            retcode += 1
            cprint("[%s]: no hints entry" % path)
            continue
        for k in ("low", "normal", "high"):
            if k not in hints:
                retcode += 1
                cprint("[hint=%s]: not present" % k)
                continue

    return retcode


#def djrepo_regmd5(options):
#    return 0


def djrepo_convert(options):
    raise NotImplementedError()
    # Version 1.0 to 2.0
    new_version = "2.0"
    for path in options.paths:
        report = DojoReport.from_file(path)
        old_version = report["version"]
        if  old_version == new_version: continue
        report = report.to_version(new_version)
        report.json_write()

    return 0


def djrepo_recalc(options):
    """
    Recompute the deltafactor from the data stored in the dojoreport.
    This function is used when the deltafactor reference results have been changed.
    """
    for path in options.paths:
        pseudo = dojopseudo_from_file(path)
        data = pseudo.dojo_report["deltafactor"]

        # Recompute delfactor for the different ecut.
        for ecut, entry in data.items():
            new_entry, eos_fit = dojo_dfact_results(pseudo, entry["num_sites"], entry["volumes"], entry["etotals"])
            print("ecut: ", ecut, "(new - old) df:", new_entry["dfact_meV"] - entry["dfact_meV"])
            data[ecut] = new_entry

        # Write new djson file.
        pseudo.dojo_report.json_write()

    return 0


def djrepo_pop(options):
    """Remove trials from djrepo files."""
    for path in options.paths:
        pseudo = dojopseudo_from_file(path)
        report = pseudo.dojo_report
        count = 0
        for trial in options.trials:
            if report.pop(trial, None) is not None:
                cprint("[%s] removed from %s" % (trial, path), "yellow")
                count += 1

        if count != 0:
            pseudo.dojo_report.json_write()


def djrepo_copyhints(options):
    """
    Copy hints from djrepo files.
    """
    if len(options.paths) != 2:
        raise ValueError("Expecting two arguments with source and destination")

    def copy_hints(from_djrepo, to_djrepo, verbose):
        from_report = DojoReport.from_file(from_djrepo)
        if "hints" not in from_report:
            print("hints are not present in source", to_djrepo)
            return 1
        to_report = DojoReport.from_file(to_djrepo)
        if "hints" in to_report:
            print("hints are already available in destination", to_djrepo)
            return 2

        if verbose:
            print("Copying hints %s --> %s" % (from_djrepo, to_djrepo))

        # Add hints, and tags. Write new json file.
        to_report["hints"] = from_report["hints"]
        to_report["tags"] = from_report.get("tags", [])
        to_report.json_write()
        return 0

    retcode = 0
    src, dest = options.paths

    if os.path.isfile(src):
        # Copy hints from file
        assert os.path.isfile(dest)
        assert src.endswith(".djrepo") and dest.endswith(".djrepo")
        return copy_hints(src, dest, options.verbose)

    # Walk through directory src and copy all hints.
    for dirpath, dirnames, filenames in os.walk(src):
        for f in filenames:
            src_path = os.path.join(dirpath, f)
            if not src_path.endswith(".djrepo"): continue
            relpath = os.path.relpath(src_path, src)
            dest_path = os.path.join(dest, relpath)
            # This to copy from SR --> FR
            dest_path = os.path.join(dest, relpath.replace(".djrepo", "_r.djrepo"))
            if not os.path.exists(dest_path):
                print("Ignoring non-existent target ", dest_path)
                continue
            rc = copy_hints(src_path, dest_path, options.verbose)
            if rc != 0:
                print("Non zero-exit status while copying hints", src_path, " --> ", dest_path)
                retcode += rc

    return retcode


@prof_main
def main():
    def str_examples():
        return """\
Example usage:
    djrepo.py check     => Check djrepo files.
    djrepo.py convert   => Convert djrepo to new format.
    djrepo.py recalc    => Recalculate deltafactor in djrepo
    djrepo.py copyhints => Copy hints from previous set of djrepo files.
"""

    def show_examples_and_exit(err_msg=None, error_code=0):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg:
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    def parse_trials(s):
        if s == "all": return ["df", "gbrv"]
        return s.split(",")

    def parse_trials(s):
        if s == "all": return DojoReport.ALL_TRIALS
        trials = s.split(",")
        unknowns = [t for t in trials if t not in DojoReport.ALL_TRIALS]
        if unknowns:
            cprint("The following names are not valid PseudoDojo trials", "red")
            print(str(unknowns))
            raise SystemExit()
        return trials

    # Parent parser for commands that need to know on which subset of pseudos we have to operate.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('paths', nargs="+", help="Pseudopotential file or directory containing pseudos")
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='verbose, can be supplied multiple times to increase verbosity')

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')

    # Subparser for check command.
    p_check = subparsers.add_parser('check', parents=[copts_parser], help=djrepo_check.__doc__)

    # Subparser for convert command.
    #p_convert = subparsers.add_parser('convert', parents=[copts_parser], help=djrepo_convert.__doc__)

    p_pop = subparsers.add_parser('pop', parents=[copts_parser], help=djrepo_pop.__doc__)
    p_pop.add_argument('--trials', default="all",  type=parse_trials, help="Select trials to remove.")

    # Subparser for fix command.
    p_fix = subparsers.add_parser('fix', parents=[copts_parser], help=djrepo_fix.__doc__)

    # Subparser for recalc command.
    p_recalc = subparsers.add_parser('recalc', parents=[copts_parser], help=djrepo_recalc.__doc__)

    # Subparser for copyhints command.
    p_copyhints = subparsers.add_parser('copyhints', parents=[copts_parser], help=djrepo_copyhints.__doc__)

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc:
        show_examples_and_exit(error_code=1)

    if options.command == "copyhints":
        return djrepo_copyhints(options)

    # Get paths
    def get_djrepo_paths(options):
        """
        Find and retur list of djson paths from options
        Accepts filepaths or directory.
        """
        paths = options.paths

        if len(paths) == 1:
            # Handle directory argument
            if os.path.isdir(paths[0]):
                top = os.path.abspath(paths[0])
                paths = find_exts(top, ["djrepo"], exclude_dirs="_*")
            # Handle glob syntax e.g. "./*.psp8"
            elif "*" in paths[0]:
                paths = [f for f in glob.glob(paths[0]) if f.endswith(".djrepo")]
        else:
            paths = [p for p in paths if p.endswith(".djrepo")]

        return paths

    options.paths = get_djrepo_paths(options)
    if not options.paths:
        cprint("Empty djrepo list. Returning", "magenta")
        return 1

    if options.verbose:
        print("Got %s djrepo files" % len(options.paths))
        if options.verbose > 1:
            for i, p in enumerate(options.paths): print("[%d] %s" % (i, p))

    # Dispatch
    return globals()["djrepo_" + options.command](options)


if __name__ == "__main__":
    sys.exit(main())

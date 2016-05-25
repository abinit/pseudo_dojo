#!/usr/bin/env python
"""Script to analyze/plot data reported in the DOJO_REPORT section."""
from __future__ import division, print_function, unicode_literals

import sys
import os
import glob
import argparse
import numpy as np

from time import gmtime, strftime
from warnings import warn
from pprint import pprint
from tabulate import tabulate
from pandas import DataFrame, concat
from monty.os.path import find_exts
from monty.termcolor import cprint 
from pymatgen.util.io_utils import ask_yesno, prompt
from pseudo_dojo.core.pseudos import dojopseudo_from_file, DojoTable
from pseudo_dojo.core.dojoreport import DojoReport
from pseudo_dojo.ppcodes.oncvpsp import OncvOutputParser


def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


def dojo_figures(options):
    """
    Create figures for a dojo table.
    Currently for all pseudos in the search space, the one with the best df per element is chosen.
    This should probably come from a dojotable eventually
    """
    pseudos = options.pseudos

    if False:
        """
        read the data from a data file instead of psp files
        """
        rows = []
        with open('data') as data_file:
            for line in data_file:
                line.rstrip('\n')
                #print(line)
                data = line.split(',')
                #print(data)
                data_dict = {'name': data[0],
                            'high_dfact_meV': float(data[1]),
                            'rell_high_dfact_meV': float(data[2]),
                            'high_dfactprime_meV': float(data[3])}
                if data[5] != 'nan':
                    data_dict['high_gbrv_bcc_a0_rel_err'] = float(data[5])
                    data_dict['high_gbrv_fcc_a0_rel_err'] = float(data[7])
                rows.append(data_dict)
    else:
	# Get data from dojoreport
	data_dojo, errors = pseudos.get_dojo_dataframe()
	if errors:
	    cprint("get_dojo_dataframe returned %s errors" % len(errors), "red")
	    if options.verbose:
		for i, e in enumerate(errors): print("[%s]" % i, e)

	# add data that is not part of the dojo report
	data_pseudo = DataFrame(columns=('nv', 'valence', 'rcmin', 'rcmax') )
	for index, p in data_dojo.iterrows():
	    outfile = p.filepath.replace('.psp8', '.out')
	    parser = OncvOutputParser(outfile)
	    parser.scan()
	    if not parser.run_completed:
		raise RuntimeError("[%s] Corrupted outfile")

	    data_pseudo.loc[index] = [parser.nv, parser.valence, parser.rc_min, parser.rc_max]
  
	data = concat([data_dojo, data_pseudo], axis=1)

    """Select entries per element"""
    grouped = data.groupby("symbol")

    rows, names = [], []
    for name, group in grouped:
    
        if False: # options.semicore
            select = group.sort("nv").iloc[-1]
        elif False: # options.valence
            select = group.sort("nv").iloc[0]
        else:
            select = group.sort("high_dfact_meV").iloc[0]

        names.append(name)

        try:
            l = {k: getattr(select, k) for k in ('name', "symbol", 'Z',
                                         'high_b0_GPa', 'high_b1', 'high_v0', 'high_dfact_meV',
                                         'high_dfactprime_meV', 'high_ecut', 'high_gbrv_bcc_a0_rel_err',
                                         'high_gbrv_fcc_a0_rel_err', 'high_ecut', 'low_phonon', 'high_phonon',
                                         'low_ecut_hint', 'normal_ecut_hint', 'high_ecut_hint',
                                         'nv', 'valence', 'rcmin', 'rcmax')}
        except AttributeError as exc:
            cprint("[%s] Exception %s" (name, exc), "magenta")
            l = {k: getattr(select, k) for k in ('name', "symbol", 'Z',
                                         'high_b0_GPa', 'high_b1', 'high_v0', 'high_dfact_meV',
                                         'high_dfactprime_meV', 'high_ecut',
                                         'low_ecut_hint', 'normal_ecut_hint', 'high_ecut_hint',
                                         'nv', 'valence', 'rcmin', 'rcmax')}
        rows.append(l)

    import matplotlib.pyplot as plt
    import matplotlib.cm as mpl_cm
    from pseudo_dojo.util.ptable_plotter import ElementDataPlotterRangefixer

    cmap = mpl_cm.cool
    color = 'black'
    cmap.set_under('w', 1.)
 
    # functions for plotting
    def rcmin(elt):
        """R_c min [Bohr]"""
        return elt['rcmin']

    def rcmax(elt):
        """R_c max [Bohr]"""
        return elt['rcmax']

    def ar(elt):
        """Atomic Radius [Bohr]"""
        return elt['atomic_radii']*0.018897161646320722

    def df(elt):
        """Delta Factor [meV / atom]"""
        return elt.get('high_dfact_meV', float('NaN'))

    def dfp(elt):
        """Delta Factor Prime"""
        return elt.get('high_dfactprime_meV', float('NaN'))

    def bcc(elt):
        """GBRV BCC [% relative error]"""
        try:
            v_bcc = elt['high_gbrv_bcc_a0_rel_err'] if str(elt['high_gbrv_bcc_a0_rel_err']) != 'nan' else -99
            # print(v_bcc)
            return v_bcc
        except KeyError:
            #print('bcc func fail: ', elt)
            return -99 #float('NaN')

    def fcc(elt):
        """GBRV FCC [% relative error]"""
        try:
            v_fcc = elt['high_gbrv_fcc_a0_rel_err'] if str(elt['high_gbrv_fcc_a0_rel_err']) != 'nan' else -99
            #print(v_fcc)
            return v_fcc 
        except KeyError:
            #print('fcc func fail: ', elt)
            return -99 #float('NaN')

    def low_phon_with(elt):
        """Acoustic mode low_cut"""
        try:
            return elt['low_phonon'][0]
        except (KeyError, TypeError):
            #print('low_phon wiht func fail: ', elt)
            return float('NaN')

    def high_phon_with(elt):
        """AC mode [\mu eV]"""
        try:
            return elt['high_phonon'][0]*1000
        except (KeyError, TypeError):
            #print('high_phon with func fail: ', elt)
            return float('NaN')
    
    def high_ecut(elt):
        """ecut high [Ha]"""
        return elt.get('high_ecut_hint', float('NaN'))

    def low_ecut(elt):
        """ecut low [Ha]"""
        return elt.get('low_ecut_hint', float('NaN'))

    def normal_ecut(elt):
        """ecut normal [Ha]"""
        return elt.get('normal_ecut_hint', float('NaN'))

    els = []
    elsgbrv = []
    elsphon = []
    rel_ers = []
    elements_data = {}

    for el in rows:
        symbol = el["symbol"]

        try:
            rel_ers.append(max(abs(el['high_gbrv_bcc_a0_rel_err']),abs(el['high_gbrv_fcc_a0_rel_err'])))
        except (TypeError, KeyError) as exc:
            if options.verbose: print(exc)

        if el['high_dfact_meV'] > 0:
            elements_data[symbol] = el
            els.append(symbol)
        else:
            cprint('[%s] failed reading high_dfact_meV %s:' % (symbol, el['high_dfact_meV']), "magenta")

        try:
            if el['high_gbrv_bcc_a0_rel_err'] > -100 and el['high_gbrv_fcc_a0_rel_err'] > -100:
                elsgbrv.append(symbol)
        except (KeyError, TypeError) as exc:
            cprint('[%s] failed reading gbrv' % symbol, "magenta")
            if options.verbose: print(exc)

        try:
            if len(el['high_phonon']) > 2:
                elsphon.append(symbol)
        except (KeyError, TypeError) as exc:
            cprint('[%s] failed reading high_phonon' % symbol, "magenta")
            if options.verbose: print(exc)
      
    try:
        max_rel_err = 0.05 * int((max(rel_ers) / 0.05) + 1)
    except ValueError:
        max_rel_err = 0.20

    # plot the GBRV/DF results periodic table
    epd = ElementDataPlotterRangefixer(elements=els, data=elements_data)
    cm1 = mpl_cm.jet
    cm2 = mpl_cm.jet
    cm1.set_under('w', 1.0)
    epd.ptable(functions=[bcc,fcc,df], font={'color': color}, cmaps=[cm1, cm1, cm2],
               #clims=[[-max_rel_err, max_rel_err],[-max_rel_err, max_rel_err], [-20,20]])
               clims=[[-0.6,0.6],[-0.6, 0.6], [-4,4]])
    plt.show()

    # Test different color maps
    #for cm2 in [mpl_cm.PiYG_r, mpl_cm.PRGn_r,mpl_cm.RdYlGn_r]:
    #     epd = ElementDataPlotterRangefixer(elements=els, data=elements_data)
    #     epd.ptable(functions=[bcc,fcc,df], font={'color':color}, cmaps=[cm1,cm1,cm2],
    #           clims=[[-max_rel_err,max_rel_err],[-max_rel_err, max_rel_err], [0,3]])
    #     plt.show()
    #plt.savefig('gbrv.eps', format='eps')

    # plot the periodic table with deltafactor and deltafactor prime.
    epd = ElementDataPlotterRangefixer(elements=els, data=elements_data)
    epd.ptable(functions=[df,dfp], font={'color':color}, cmaps=cmap, clims=[[0, 6]])
    plt.show()
    #plt.savefig('df.eps', format='eps')

    # plot the GBVR results periodic table
    epd = ElementDataPlotterRangefixer(elements=elsgbrv, data=elements_data)
    epd.ptable(functions=[bcc,fcc], font={'color':color}, cmaps=mpl_cm.jet, clims=[[-max_rel_err, max_rel_err]])
    plt.show()
    #plt.savefig('gbrv.eps', format='eps')

    # plot the hints periodic table
    epd = ElementDataPlotterRangefixer(elements=els, data=elements_data)
    cm = mpl_cm.cool
    cm.set_under('w', 1.0)
    epd.ptable(functions=[low_ecut, high_ecut, normal_ecut], font={'color':color}, clims=[[6, 80]],  cmaps=cmap)
    plt.show()
    #plt.savefig('rc.eps', format='eps')

    # plot the radii periodic table
    epd = ElementDataPlotterRangefixer(elements=els, data=elements_data)
    epd.ptable(functions=[rcmin, rcmax, ar], font={'color':color}, clims=[[0, 4]], cmaps=cmap)
    plt.show()
    #plt.savefig('rc.eps', format='eps')

    # plot the accoustic mode periodic table
    epd = ElementDataPlotterRangefixer(elements=elsphon, data=data)
    cm = mpl_cm.winter
    cm.set_under('orange', 1.0)
    epd.ptable(functions=[high_phon_with], font={'color':color}, cmaps=cm, clims=[[-2, 0]])
    plt.show()
    #plt.savefig('rc.eps', format='eps')

    return 0


def dojo_plot(options):
    """Plot DOJO results for specified pseudos."""
    pseudos = options.pseudos

    for pseudo in pseudos:
        if not pseudo.has_dojo_report:
            cprint("%s does not contain the DOJO_REPORT section" % pseudo.filepath, "magenta")
            continue

        report = pseudo.dojo_report

        if options.verbose:
            print(pseudo)
	    try:
		error = report.check()
		if error:
		    cprint("[%s] Validation error" % pseudo.basename, "red")
                    print(error)
                    print("")
	    except Exception as exc:
		cprint("[%s] Python exception in report_check" % pseudo.basename, "red")
		print(exc)

        # ebands
        if report.has_trial("ebands") and any(k in options.what_plot for k in ("all", "ebands")):
            try:
                report.plot_ebands(title=pseudo.basename)
            except Exception as exc:
                cprint(exc, "red")

        # Deltafactor
        if report.has_trial("deltafactor") and any(k in options.what_plot for k in ("all", "df")):
            if options.eos: 
                # Plot EOS curve
                report.plot_deltafactor_eos(title=pseudo.basename)

            # Plot total energy convergence.
            fig = report.plot_deltafactor_convergence(xc=pseudo.xc, title=pseudo.basename)
            fig = report.plot_etotal_vs_ecut(title=pseudo.basename)

        # GBRV
        if any(k in options.what_plot for k in ("all", "gbrv")):
            count = 0
            for struct_type in ("fcc", "bcc"):
                trial = "gbrv_" + struct_type
                if report.has_trial(trial):
                    count += 1
                    if options.eos:
                        report.plot_gbrv_eos(struct_type=struct_type, title=pseudo.basename)
            if count:
                report.plot_gbrv_convergence(title=pseudo.basename)

        # phonon
        if any(k in options.what_plot for k in ("all", "phonon")):
            if report.has_trial("phonon"):
                report.plot_phonon_convergence(title=pseudo.basename)

        #if any(k in options.what_plot for k in ("all", "phwoa")):
        #    if report.has_trial("phwoa"):
        #        report.plot_phonon_convergence(title=pseudo.basename+"W/O ASR", woasr=True)

    return 0


def dojo_notebook(options):
    """
    Generate an ipython notebook for each pseudopotential and open it in the browser.
    """
    from pseudo_dojo.util.notebook import make_open_notebook
    retcode = 0
    for p in options.pseudos:
        retcode += make_open_notebook(p.filepath)
        if retcode != 0: break

    return retcode


def dojo_compare(options):
    """Compare DOJO results for multiple pseudos."""
    pseudos = options.pseudos
    for z in pseudos.zlist: 
        pseudos_z = pseudos[z]
        if len(pseudos_z) > 1:
            pseudos_z.dojo_compare(what=options.what_plot)
        else:
            print("Found only one pseudo for Z=%s" % z)

    return 0


def dojo_trials(options):
    """Visualize the results of the different tests."""
    pseudos = options.pseudos

    # Build pandas DataFrame
    data, errors = pseudos.get_dojo_dataframe()
    #print(data)

    if errors:
        cprint("ERRORS:", "red")
        pprint(errors)
    
    #import matplotlib.pyplot as plt
    #data.plot_trials(savefig=options.savefig)
    #data.plot_hist(savefig=options.savefig)
    #data.sns_plot(savefig=options.savefig)
    #print(data["high_dfact_meV"])
    import matplotlib.pyplot as plt
    ax = data["high_dfact_meV"].hist()
    ax.set_xlabel("Deltafactor [meV]")
    plt.show()
    return 0


def dojo_table(options):
    """Build and show a pandas table."""
    pseudos = options.pseudos
    data, errors = pseudos.get_dojo_dataframe()

    if errors:
        cprint("get_dojo_dataframe returned %s errors" % len(errors), "red")
        if options.verbose:
            for i, e in enumerate(errors): print("[%s]" % i, e)

    #data.tabulate()
    #print(data.columns)
    #frame = data[["high_dfact_meV", "high_dfactprime_meV", "high_gbrv_bcc_a0_rel_err", "high_gbrv_fcc_a0_rel_err"]]
    #data["SOC"] = ["_r" in s for s in data["filepath"]]
    #data = data[["SOC", "high_dfact_meV", "high_gbrv_bcc_a0_rel_err", "high_gbrv_fcc_a0_rel_err", "Z"]]
    #print(tabulate(data, headers="keys"))
    #import matplotlib.pyplot as plt
    #import seaborn as sns
    #g = sns.FacetGrid(data, hue="SOC")
    #g.map(plt.scatter, "Z", "high_dfact_meV"); g.add_legend(); plt.show()
    #g.map(plt.scatter, "Z", "high_gbrv_bcc_a0_rel_err"); g.add_legend(); plt.show()
    #g.map(plt.scatter, "Z", "high_gbrv_fcc_a0_rel_err"); g.add_legend(); plt.show()
    #return

    if False:
        """Select best entries"""
        #best = {}
        #symbols = set(data["symbol"])
        #for sym in symbols:
        grouped = data.groupby("symbol")
        #print(grouped.groups)

        rows, names = [], []
        for name, group in grouped:
            #print(name, group["high_dfact_meV"])
            best = group.sort("high_dfact_meV").iloc[0]
            names.append(name)
            #print(best.keys())
            #print(best.name)
            l = {k: getattr(best, k) for k in ("name", "Z", 'high_b0_GPa', 'high_b1', 'high_v0', 
                                               'high_dfact_meV', 'high_ecut_deltafactor')} 
            rows.append(l)

        import pandas
        best_frame = pandas.DataFrame(rows, index=names)
        best_frame = best_frame.sort("Z")
        print(tabulate(best_frame, headers="keys"))
        print(tabulate(best_frame.describe(),  headers="keys"))

        import matplotlib.pyplot as plt
        best_frame["high_dfact_meV"].hist(bins=100)
        plt.show()
        return 0

    if errors:
        cprint("ERRORS:", "red")
        pprint(errors)

    #accuracies = ["normal", "high"]
    accuracies = ["low", "normal", "high"]
    keys = ["dfact_meV", "dfactprime_meV", "v0", "b0_GPa", "b1", "ecut_deltafactor", "ecut_hint"]
    columns = ["symbol"] + [acc + "_" + k for k in keys for acc in accuracies]
    #print(columns)

    #data = data[data["high_dfact_meV"] <= data["high_dfact_meV"].mean()]
    #data = data[data["high_dfact_meV"] <= 9]

    try:
        data["low_dfact_abserr"] = data["low_dfact_meV"] - data["high_dfact_meV"]
        data["normal_dfact_abserr"] =  data["normal_dfact_meV"] - data["high_dfact_meV"]
        data["low_dfact_rerr"] = 100 * (data["low_dfact_meV"] - data["high_dfact_meV"]) / data["high_dfact_meV"]
        data["normal_dfact_rerr"] = 100 * (data["normal_dfact_meV"] - data["high_dfact_meV"]) / data["high_dfact_meV"]

        for k in ["v0", "b0_GPa", "b1"]:
            data["low_" + k + "_abserr"] = data["low_" + k] - data["high_" + k]
            data["normal_" + k + "_abserr"] = data["normal_" + k] - data["high_" + k]
            data["low_" + k + "_rerr"] = 100 * (data["low_" + k] - data["high_" + k]) / data["high_" + k]
            data["normal_" + k + "_rerr"] = 100 * (data["normal_" + k] - data["high_" + k]) / data["high_" + k]
    except Exception as exc:
        cprint("Python exception: %s" % type(exc), "red")
        if options.verbose: print(exc)
 
    try:
        for acc in ['low', 'normal', 'high']:
            data[acc + "_abs_fcc"] = abs(data[acc + "_gbrv_fcc_a0_rel_err"])
            data[acc + "_abs_bcc"] = abs(data[acc + "_gbrv_bcc_a0_rel_err"])
    except KeyError:
        cprint('no GBRV data', "magenta")

    wrong = data[data["high_b1"] < 0]
    if not wrong.empty:
        cprint("WRONG".center(80, "*"), "red")
        print(wrong)

    #data = calc_errors(data)
    #data.to_json('table.json')

    try:
        data = data[
                 [acc + "_dfact_meV" for acc in accuracies]
               + [acc + "_dfactprime_meV" for acc in accuracies]
               + [acc + "_ecut_deltafactor" for acc in accuracies]
               + [acc + "_gbrv_fcc_a0_rel_err" for acc in accuracies]
               + [acc + "_gbrv_bcc_a0_rel_err" for acc in accuracies]
               + [acc + "_abs_fcc" for acc in accuracies]
               + [acc + "_abs_bcc" for acc in accuracies]
               + [acc + "_ecut_hint" for acc in accuracies]
                   ]
    except KeyError:
        data = data[
                 [acc + "_dfact_meV" for acc in accuracies]
               + [acc + "_ecut_deltafactor" for acc in accuracies]
               + [acc + "_dfactprime_meV" for acc in accuracies]
               + [acc + "_ecut_hint" for acc in accuracies]
                   ]

    print("\nONCVPSP TABLE:\n") #.center(80, "="))
    tablefmt = "grid"
    floatfmt = ".3f"

    accuracies = ['low', 'high']
    columns = [acc + "_dfact_meV" for acc in accuracies] 
    columns += [acc + "_ecut_deltafactor" for acc in accuracies] 

    print(tabulate(data[columns], headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))
    if len(data) > 5: 
        print(tabulate(data[columns].describe(), headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))

    accuracies = ['low', 'high']
    columns = [acc + "_dfactprime_meV" for acc in accuracies]
    columns += [acc + "_ecut_deltafactor" for acc in accuracies]

    print(tabulate(data[columns], headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))
    if len(data) > 5:
        print(tabulate(data[columns].describe(), headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))

    try:
        columns = [acc + "_gbrv_fcc_a0_rel_err" for acc in accuracies]
        columns += [acc + "_gbrv_bcc_a0_rel_err" for acc in accuracies]
        columns += [acc + "_abs_fcc" for acc in accuracies]
        columns += [acc + "_abs_bcc" for acc in accuracies]
        print(tabulate(data[columns], headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))
        if len(data) > 5:
            print(tabulate(data[columns].describe(), headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))
    except KeyError as exc:
        cprint('No GBRV data', "red")    
        if options.verbose: print("Python exception:\n", str(exc))

    accuracies = ['low', 'normal', 'high']
    columns = [acc + "_ecut_hint" for acc in accuracies]
    print(tabulate(data[columns], headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))
    if len(data) > 5:
        print(tabulate(data[columns].describe(), headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))

    #print(data.to_string(columns=columns))

    if len(data) > 5: 
        bad = data[data["high_dfact_meV"] > data["high_dfact_meV"].mean() + data["high_dfact_meV"].std()]
        good = data[data["high_dfact_meV"] < data["high_dfact_meV"].mean()]
        print("\nPSEUDOS with high_dfact > mean plus one std (%.3f + %.3f):\n" % (
              data["high_dfact_meV"].mean(), data["high_dfact_meV"].std())) # ".center(80, "*"))
        print(tabulate(bad[["high_dfact_meV", "high_ecut_deltafactor"]], headers="keys", 
             tablefmt=tablefmt, floatfmt=floatfmt))

    #gbrv_fcc_bad = data[data["high_gbrv_fcc_a0_rerr"] > (data["high_gbrv_fcc_a0_rerr"].abs()).mean()]
    #print("\nPSEUDOS with high_dfact > mean:\n") # ".center(80, "*"))
    #print(tabulate(bad, headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))

    return 0


def dojo_dist(options):
    """
    Plot the distribution of the deltafactor and of the relative error for the GBRV fcc/bcc tests.
    """
    fig = options.pseudos.plot_dfgbrv_dist()
    return 0


def dojo_check(options):
    """Check validity of pseudodojo report."""
    retcode = 0
    for p in options.pseudos:
        try:
            report = p.dojo_report
        except Exception as exc:
            cprint("[%s] Invalid dojo_report" % p.basename, "red")
            if options.verbose:
                print("Python Exception:\n%s", exc)
            retcode += 1
            continue

        #print(report)
        # Comment this to fix the md5 checksum in the pseudos
        #p.check_and_fix_dojo_md5()

        #if "ppgen_hints" not in report: # and "deltafactor" not in report:
        #    print(p.basename, "old version without ppgen_hints")
        #    continue

        try:
            error = report.check(check_trials=options.check_trials)
            if error:
                retcode += 1
                if options.verbose:
		    cprint("[%s] Validation error" % p.basename, "red")
                    print(error)
                    print("")
		else:
		    cprint("[%s] Validation error. Use -v for more info" % p.basename, "red")

        except Exception as exc:
            retcode += 1
            cprint("[%s] Python exception:" % p.basename, "red")
            if options.verbose:
                print(str(exc))

    return retcode


def dojo_make_hints(options):
    """Add hints for energy cutoffs"""
    for pseudo in options.pseudos:

        if not pseudo.has_dojo_report:
            cprint("[%s] No DojoReport. Ignoring pseudo" % p.basename, "red")
            continue

        report = pseudo.dojo_report

        try:
            hints = report.compute_hints()
            print("hints for %s computed from deltafactor prime: %s" % (pseudo.basename, hints))
            report.plot_deltafactor_convergence(xc=pseudo.xc)
 
            ans = ask_yesno("Do you accept the hints? [Y]")
            if ans:
                report.add_hints(hints)
                print(report.has_hints)
                report.json_write(pseudo.djrepo_path)
            else:
                print("The dojoreport contains ecuts :\n%s" % report.ecuts)
                new_ecuts = prompt("Enter new ecuts to compute (comma-separated values or empty string to abort)")
                if not new_ecuts:
                    print("Exit requested by user")
                    return 
                new_ecuts = np.array([float(k) for k in new_ecuts.strip().split(",")])
                print("new_ecuts", new_ecuts)
                report.add_ecuts(new_ecuts)
                report.json_write(pseudo.djrepo_path)

        except Exception as exc:
            cprint("[%s]: python exception: %s" % (pseudo.basename, type(exc)), "red")
            if options.verbose: print(straceback())

    return 0


def dojo_validate(options):
    """Validate the pseudo."""
    pseudos = options.pseudos
    data, errors = pseudos.get_dojo_dataframe()

    for p in pseudos:
        try:
            # test if report is present
            if not p.has_dojo_report:
                cprint("[%s] No DojoReport. Ignoring pseudo" % p.basename, "red")
                continue
            report = p.dojo_report

            # test if already validated
            if 'validation' in report:
                cprint('this pseudo was validated by %s on %s.' % (
                       report['validation']['validated_by'], report['validation']['validated_on']), "red")
		if not ask_yesno('Would you like to validate it again? [Y/n]'):
		    continue

            # test for ghosts
            print('\n= GHOSTS TEST ===========================================\n')
            if report.has_trial('ebands'):
                for ecut in report['ebands']:
                    if "ghost_free_upto_eV" in report["ebands"][ecut]:
                        print('%s: Pseudo is reported to be ghost free up to %s eV' % (
                          ecut, report["ebands"][ecut]["ghost_free_upto_eV"]))
                    else:
                        report.plot_ebands(ecut=ecut)
                        ans = float(prompt('Please enter the energy (eV) up until there is no sign of ghosts:\n'))
                        if ans > 0:
                            report["ebands"][ecut]["ghost_free_upto_eV"] = ans
                            report.json_write(p.djrepo_path)
            else:
                cprint('no ebands trial present, pseudo cannot be validated', "red")
                continue

            # test trials
            print('\n= TRIALS TEST ===========================================\n')
            try:
                error = report.check()
                if error:
                    cprint("[%s] Validation problem" % p.basename, "red")
                    if options.verbose:
                        print(error)
                        print()

            except Exception as exc:
                cprint("[%s] Python exception: %s" % (p.basename, type(exc)), "red")
                if options.verbose:
                    print(exc)
                    print("")

            accuracies = ["low", "normal", "high"]
            keys = ['hint', 'deltafactor', 'gbrv_bcc', 'gbrv_fcc', 'phonon']

            tablefmt = "grid"
            floatfmt=".2f"

            print('\n= ECUTS  TEST ===========================================\n')
            for acc in accuracies:
                columns = ["symbol"]
                headers = ["symbol"]
                for k in keys:
                    entry = acc + "_ecut_" + k 
                    if entry in data:
                        columns.append(entry)
                        headers.append(k)
                print('ECUTS for accuracy %s:' % acc)
                print(tabulate(data[columns], headers=headers, tablefmt=tablefmt, floatfmt=floatfmt))

            # test hash
            # plot the model core charge
           
        except Exception as exc:
            cprint("[%s] python exception" % p.basename, "red")
            if options.verbose: print(straceback())
        
        # ask the final question
        if ask_yesno('Will you validate this pseudo? [n]'):
            name = prompt("Please enter your name for later reference: ")
            report['validation'] = {'validated_by': name, 'validated_on': strftime("%Y-%m-%d %H:%M:%S", gmtime())}
            report.json_write(p.djrepo_path)

    return 0


def main():
    def str_examples():
        return """\
Usage example:
    dojodata.py plot H.psp8                ==> Plot dojo data for pseudo H.psp8
    dojodata.py compare H.psp8 H-low.psp8  ==> Plot and compare dojo data for pseudos H.psp8 and H-low.psp8
    dojodata.py trials H.psp8 -r 1
    dojodata.py table .                    ==> Build table (find all psp8 files within current directory)
    dojodata.py figures .                  ==> Plot periodic table figures
    dojodata.py notebook H.psp8            ==> Generate ipython notebook and open it in the browser
    dojodata.py check table/*/*_r.psp8 -v --check-trials=gbrv_fcc,gbrv_bcc
"""

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    def parse_rows(s):
        if not s: return []
        tokens = s.split(",")
        return list(map(int, tokens)) if tokens else []

    def parse_symbols(s):
        if not s: return []
        return s.split(",")

    # Parent parser for commands that need to know on which subset of pseudos we have to operate.
    pseudos_selector_parser = argparse.ArgumentParser(add_help=False)
    pseudos_selector_parser.add_argument('pseudos', nargs="+", help="Pseudopotential file or directory containing pseudos")
    pseudos_selector_parser.add_argument('-s', "--symbols", type=parse_symbols, 
        help=("List of chemical symbols to include or exclude."
              "Example --symbols=He,Li to include He and Li, --symbols=-He to exclude He"))
    pseudos_selector_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='Verbose, can be supplied multiple times to increase verbosity')

    # Options for pseudo selection.
    group = pseudos_selector_parser.add_mutually_exclusive_group()
    group.add_argument("-r", '--rows', default="", type=parse_rows, help="Select these rows of the periodic table.")
    group.add_argument("-f", '--family', type=str, default="", help="Select this family of the periodic table.")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")
    parser.add_argument('--seaborn', action="store_true", help="Use seaborn settings")

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    plot_options_parser = argparse.ArgumentParser(add_help=False)
    plot_options_parser.add_argument("-w", "--what-plot", type=str, default="all", 
                                      help="Quantity to plot e.g df for deltafactor, gbrv for GBRV tests")
    plot_options_parser.add_argument("-e", "--eos", action="store_true", help="Plot EOS curve")

    # Subparser for plot command.
    p_plot = subparsers.add_parser('plot', parents=[pseudos_selector_parser, plot_options_parser],
                                   help=dojo_plot.__doc__)

    # Subparser for notebook command.
    p_notebook = subparsers.add_parser('notebook', parents=[pseudos_selector_parser],
                                       help=dojo_notebook.__doc__)

    # Subparser for compare.
    p_compare = subparsers.add_parser('compare', parents=[pseudos_selector_parser, plot_options_parser],
                                      help=dojo_compare.__doc__)

    # Subparser for figures
    p_figures = subparsers.add_parser('figures', parents=[pseudos_selector_parser], help=dojo_figures.__doc__)

    # Subparser for table command.
    p_table = subparsers.add_parser('table', parents=[pseudos_selector_parser], help=dojo_table.__doc__)

    # Subparser for dist command.
    p_dist = subparsers.add_parser('dist', parents=[pseudos_selector_parser], help=dojo_dist.__doc__)

    # Subparser for trials command.
    p_trials = subparsers.add_parser('trials', parents=[pseudos_selector_parser], help=dojo_trials.__doc__)
    p_trials.add_argument("--savefig", type=str, default="", help="Save plot to savefig file")

    # Subparser for check command.
    def parse_trials(s):
        if s == "all": return DojoReport.ALL_TRIALS
        return s.split(",")

    p_check = subparsers.add_parser('check', parents=[pseudos_selector_parser], help=dojo_check.__doc__)
    p_check.add_argument("--check-trials", type=parse_trials, default="all", help="List of trials to check")

    # Subparser for validate command.
    p_validate = subparsers.add_parser('validate', parents=[pseudos_selector_parser], help=dojo_validate.__doc__)

    # Subparser for make_hints command.
    p_make_hints = subparsers.add_parser('make_hints', parents=[pseudos_selector_parser],
                                         help=dojo_make_hints.__doc__)

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc: 
        show_examples_and_exit(error_code=1)

    # loglevel is bound to the string value obtained from the command line argument. 
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    def get_pseudos(options):
        """
        Find pseudos in paths, return :class:`DojoTable` object sorted by atomic number Z.
        Accepts filepaths or directory.
        """
        exts = ("psp8",)

        paths = options.pseudos

        if len(paths) == 1:
            # Handle directory argument
            if os.path.isdir(paths[0]):
                top = os.path.abspath(paths[0])
                paths = find_exts(top, exts, exclude_dirs="_*")
            # Handle glob syntax e.g. "./*.psp8" 
            elif "*" in paths[0]:
                paths = glob.glob(paths[0])

        if options.verbose > 1: print("Will read pseudo from: %s" % paths)

        pseudos = []
        for p in paths:
            try:
                pseudo = dojopseudo_from_file(p)
                if pseudo is None: 
                    cprint("[%s] Pseudo.from_file returned None. Something wrong in file!" % p, "red")
                    continue
                pseudos.append(pseudo)

            except Exception as exc:
                cprint("[%s] Python exception. This pseudo will be ignored" % p, "red")
                if options.verbose: print(exc)

        table = DojoTable(pseudos)

        # Here we select a subset of pseudos according to family or rows
        if options.rows:
            table = table.select_rows(options.rows)
        elif options.family:
            table = table.select_families(options.family)

        # here we select chemical symbols.
        if options.symbols:
            table = table.select_symbols(options.symbols)

        return table.sort_by_z()

    # Build DojoTable from the paths specified by the user.
    options.pseudos = get_pseudos(options)
    if not options.pseudos:
	cprint("Empty pseudopotential list. Returning", "magenta")
	return 1

    if options.seaborn:
        import seaborn as sns
        #sns.set(style='ticks', palette='Set2')
        sns.set(style="dark", palette="Set2")
        #And to remove "chartjunk", do:
        #sns.despine()
        #plt.tight_layout()
        #sns.despine(offset=10, trim=True)

    # Dispatch
    return globals()["dojo_" + options.command](options)


if __name__ == "__main__":
    try:
        do_prof = sys.argv[1] == "prof"
        if do_prof: sys.argv.pop(1)
    except Exception: 
        do_prof = False

    if do_prof:
        import pstats, cProfile
        cProfile.runctx("main()", globals(), locals(), "Profile.prof")
        s = pstats.Stats("Profile.prof")
        s.strip_dirs().sort_stats("time").print_stats()
    else:
        sys.exit(main())

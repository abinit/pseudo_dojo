#!/usr/bin/env python
"""Script to analyze/plot data reported in the DOJO_REPORT section."""
import sys
import os
import glob
import argparse
import numpy as np

from pprint import pprint
from tabulate import tabulate
#from pandas import DataFrame, concat
from monty.os.path import find_exts
from monty.functools import prof_main
from monty import termcolor
from monty.termcolor import cprint
from pseudo_dojo.core.pseudos import dojopseudo_from_file, DojoTable
#from pseudo_dojo.ppcodes.oncvpsp import OncvOutputParser
from pseudo_dojo.pseudos import check_pseudo


def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


def prompt(question):
    my_input = input
    return my_input(question)


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
        print('test')

        if errors:
            cprint("get_dojo_dataframe returned %s errors" % len(errors), "red")
            if not options.verbose:
                print("Use --verbose for details.")
            else:
                for i, e in enumerate(errors):
                    print("[%s]" % i, e)

        # add data that is not part of the dojo report
#        data_pseudo = DataFrame(columns=('nv', 'valence', 'rcmin', 'rcmax') )
#        for index, p in data_dojo.iterrows():
#            print(p.keys())
#            outfile = p.filepath.replace('.psp8', '.out')
#            parser = OncvOutputParser(outfile)
#            parser.scan()
#            if not parser.run_completed:
#                raise RuntimeError("[%s] Corrupted outfile")

#            data_pseudo.loc[index] = [parser.nv, parser.valence, parser.rc_min, parser.rc_max]

#        data = concat([data_dojo, data_pseudo], axis=1)
        data = data_dojo

    # Select "best" entries per element.
    rows, names = [], []
    sortby, ascending = "high_dfact_meV", True

    for name, group in data.groupby("symbol"):
        # Sort group and select best pseudo depending on sortby and ascending.
        select = group.sort_values(sortby, ascending=ascending).iloc[0]
        l = {k: getattr(select, k, None) for k in (
                                             'name', "symbol", 'Z',
                                             'high_b0_GPa', 'high_b1', 'high_v0', 'high_dfact_meV',
                                             'high_dfactprime_meV', 'high_ecut', 'high_gbrv_bcc_a0_rel_err',
                                             'high_gbrv_fcc_a0_rel_err', 'high_ecut',
                                             'low_phonon', 'normal_phonon', 'high_phonon',
                                             'low_ecut_hint', 'normal_ecut_hint', 'high_ecut_hint',
                                             'nv', 'valence', 'rcmin', 'rcmax')}
        for k, v in l.items():
            if v is None: cprint("[%s] Got None for %s" % (name, k), "red")

        names.append(name)
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
        return elt['atomic_radii'] * 0.018897161646320722

    def df(elt):
        """Delta Factor [meV / atom]"""
        return elt.get('high_dfact_meV', float('NaN'))

    def dfp(elt):
        """Delta Factor Prime"""
        return elt.get('high_dfactprime_meV', float('NaN'))

    def bcc(elt):
        """GBRV BCC [% relative error]"""
        try:
            return elt['high_gbrv_bcc_a0_rel_err'] if str(elt['high_gbrv_bcc_a0_rel_err']) != 'nan' else -99
        except KeyError:
            #print('bcc func fail: ', elt)
            return float('NaN')

    def fcc(elt):
        """GBRV FCC [% relative error]"""
        try:
            return elt['high_gbrv_fcc_a0_rel_err'] if str(elt['high_gbrv_fcc_a0_rel_err']) != 'nan' else -99
        except KeyError:
            #print('fcc func fail: ', elt)
            return float('NaN')

    def low_phon_with(elt):
        """Acoustic mode low_cut"""
        try:
            return elt['low_phonon'][0]
        except (KeyError, TypeError):
            #print('low_phon wiht func fail: ', elt)
            return float('NaN')

    def high_phon_with(elt):
        r"""AC mode [\mu eV]"""
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
    #elsphon = []
    rel_ers = []
    elements_data = {}

    for el in rows:
        symbol = el["symbol"]

        # Prepare data for deltafactor
        if el['high_dfact_meV'] is None:
            cprint('[%s] failed reading high_dfact_meV %s:' % (symbol, el['high_dfact_meV']), "magenta")
        else:
            if el['high_dfact_meV'] < 0:
                cprint('[%s] negative high_dfact_meV %s:' % (symbol, el['high_dfact_meV']), "red")
                print(symbol, el['high_dfact_meV'])
            #assert el['high_dfact_meV'] >= 0
            elements_data[symbol] = el
            els.append(symbol)

        # Prepare data for GBRV
        try:
            rel_ers.append(max(abs(el['high_gbrv_bcc_a0_rel_err']), abs(el['high_gbrv_fcc_a0_rel_err'])))
        except (TypeError, KeyError) as exc:
            cprint('[%s] failed reading high_gbrv:' % symbol, "magenta")
            if options.verbose: print(exc)

        try:
            if el['high_gbrv_bcc_a0_rel_err'] > -100 and el['high_gbrv_fcc_a0_rel_err'] > -100:
                elsgbrv.append(symbol)
        except (KeyError, TypeError) as exc:
            cprint('[%s] failed reading GBRV data for ' % symbol, "magenta")
            if options.verbose: print(exc)

        #try:
        #    if len(el['high_phonon']) > 2:
        #        elsphon.append(symbol)
        #except (KeyError, TypeError) as exc:
        #    cprint('[%s] failed reading high_phonon' % symbol, "magenta")
        #    if options.verbose: print(exc)

        #if symbol == "Br":
        #    print (elements_data[symbol])

    try:
        max_rel_err = 0.05 * int((max(rel_ers) / 0.05) + 1)
    except ValueError:
        max_rel_err = 0.20

    # plot the GBRV/DF results periodic table
    epd = ElementDataPlotterRangefixer(elements=els, data=elements_data)
    cm1 = mpl_cm.jet
    cm2 = mpl_cm.jet
    cm1.set_under('w', 1.0)
    epd.ptable(functions=[bcc, fcc, df], font={'color': color}, cmaps=[cm1, cm1, cm2],
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
    epd.ptable(functions=[df, dfp], font={'color': color}, cmaps=cmap, clims=[[0, 6]])
    plt.show()
    #plt.savefig('df.eps', format='eps')

    # plot the GBVR results periodic table
    epd = ElementDataPlotterRangefixer(elements=elsgbrv, data=elements_data)
    epd.ptable(functions=[bcc, fcc], font={'color': color}, cmaps=mpl_cm.jet, clims=[[-max_rel_err, max_rel_err]])
    plt.show()
    #plt.savefig('gbrv.eps', format='eps')

    # plot the hints periodic table
    epd = ElementDataPlotterRangefixer(elements=els, data=elements_data)
    cm = mpl_cm.cool
    cm.set_under('w', 1.0)
    epd.ptable(functions=[low_ecut, high_ecut, normal_ecut], font={'color': color}, clims=[[6, 80]],  cmaps=cmap)
    plt.show()
    #plt.savefig('rc.eps', format='eps')

    # plot the radii periodic table
    epd = ElementDataPlotterRangefixer(elements=els, data=elements_data)
    epd.ptable(functions=[rcmin, rcmax, ar], font={'color': color}, clims=[[0, 4]], cmaps=cmap)
    plt.show()
    #plt.savefig('rc.eps', format='eps')

    # plot the acoustic mode periodic table
    #epd = ElementDataPlotterRangefixer(elements=elsphon, data=data)
    #cm = mpl_cm.winter
    #cm.set_under('orange', 1.0)
    #epd.ptable(functions=[high_phon_with], font={'color':color}, cmaps=cm, clims=[[-2, 0]])
    #plt.show()
    #plt.savefig('rc.eps', format='eps')

    return 0


def dojo_plot(options):
    """Plot DOJO results for specified pseudos."""
    pseudos = options.pseudos
    socs = [False, True]

    for pseudo in pseudos:
        if not pseudo.has_dojo_report:
            cprint("%s does not contain the DOJO_REPORT section" % pseudo.filepath, "magenta")
            continue

        report = pseudo.dojo_report

        if options.verbose:
            print(pseudo)
            try:
                error = report.check(check_trias=None)
                if error:
                    cprint("[%s] Validation error" % pseudo.basename, "red")
                    print(error)
                    print("")
            except Exception as exc:
                cprint("[%s] Python exception in report_check" % pseudo.basename, "red")
                print(exc)

        # ghosts
        if any(k in options.what_plot for k in ("all", "ghosts")):
            for with_soc in socs:
                report.plot_ebands(with_soc=with_soc, title=pseudo.basename)

        # Deltafactor
        if any(k in options.what_plot for k in ("all", "df")):
            if options.eos:
                # Plot EOS curve
                for with_soc in socs:
                    report.plot_deltafactor_eos(with_soc=with_soc, title=pseudo.basename)

            # Plot total energy convergence.
            for with_soc in socs:
                report.plot_deltafactor_convergence(pseudo.xc, with_soc=with_soc, title=pseudo.basename)

            for with_soc in socs:
                report.plot_etotal_vs_ecut(with_soc=with_soc, title=pseudo.basename)

        # GBRV
        if any(k in options.what_plot for k in ("all", "gbrv")):
            #print("trials", report.trials)
            count = 0
            for struct_type in ("fcc", "bcc"):
                trial = "gbrv_" + struct_type
                trial_soc = trial + "_soc"
                if report.has_trial(trial) or report.has_trial(trial_soc):
                    count += 1
                    if options.eos:
                        report.plot_gbrv_eos(struct_type=struct_type, title=pseudo.basename)
            if count:
                for with_soc in socs:
                    report.plot_gbrv_convergence(with_soc=with_soc, title=pseudo.basename)

        # phgamma
        if any(k in options.what_plot for k in ("all", "phgamma")):
            for with_soc in socs:
                report.plot_phonon_convergence(with_soc=with_soc, title=pseudo.basename)

    return 0


def dojo_nbtable(options):
    """
    Generate an ipython notebook for a pseudopotential table and open it in the browser.
    """
    return options.pseudos.make_open_notebook()


def dojo_notebook(options):
    """
    Generate an ipython notebook for each pseudopotential and open it in the browser.
    """
    from pseudo_dojo.util.notebook import make_open_notebook
    for p in options.pseudos:
        make_open_notebook(p.filepath, with_validation=not options.no_validation,
                           with_eos=True, hide_code=options.hide_code, tmpfile=not options.no_tmp)
    return 0


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


def dojo_nbcompare(options):
    """Generate ipython notebooks to compare DOJO results for multiple pseudos."""
    pseudos = options.pseudos
    for z in pseudos.zlist:
        pseudos_z = pseudos[z]
        if len(pseudos_z) > 1:
            pseudos_z.dojo_nbcompare(what=options.what_plot)
        else:
            print("Found only one pseudo for Z=%s" % z)

        return 0


def dojo_trials(options):
    """Visualize the results of the different tests."""
    pseudos = options.pseudos

    # Build pandas DataFrame
    data, errors = pseudos.get_dojo_dataframe()
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

    # Compare FR with SR pseudos.
    #pseudos.plot_scalar_vs_fully_relativistic(what="df")
    #pseudos.plot_scalar_vs_fully_relativistic(what="gbrv")
    #return 0

    #data.plot_hist()
    #data.plot_trials()
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
    #pseudos.plot_hints()
    #return 0

    if options.best:
        print("Selecting best pseudos according to deltafactor")
        best_frame = data.select_best()
        if options.json: best_frame.to_json('table.json')

        print(tabulate(best_frame, headers="keys"))
        print(tabulate(best_frame.describe(), headers="keys"))
        #best_frame["high_dfact_meV"].hist(bins=100)
        #import matplotlib.pyplot as plt
        #plt.show()
        return 0

    accuracies = ["normal", "high"]
    keys = ["dfact_meV", "dfactprime_meV", "v0", "b0_GPa", "b1", "ecut_deltafactor", "ecut_hint"]
    if options.json:
        accuracies = ["low", "normal", "high"]
        keys.append("phonon")

    columns = ["symbol"] + [acc + "_" + k for k in keys for acc in accuracies]

    #data = data[data["high_dfact_meV"] <= data["high_dfact_meV"].mean()]
    #data = data[data["high_dfact_meV"] <= 9]

    try:
        data["low_dfact_abserr"] = data["low_dfact_meV"] - data["high_dfact_meV"]
        data["normal_dfact_abserr"] = data["normal_dfact_meV"] - data["high_dfact_meV"]
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
        for acc in accuracies:
            data[acc + "_abs_fcc"] = abs(data[acc + "_gbrv_fcc_a0_rel_err"])
            data[acc + "_abs_bcc"] = abs(data[acc + "_gbrv_bcc_a0_rel_err"])
    except KeyError:
        cprint('no GBRV data', "magenta")

    try:
        wrong = data[data["high_b1"] < 0]
        if not wrong.empty:
            cprint("WRONG".center(80, "*"), "red")
            print(wrong)
    except Exception as exc:
        print(exc)

    if options.json: data.to_json('table.json')

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

    print("\nONCVPSP TABLE:\n")
    tablefmt = "grid"
    floatfmt = ".2f"

    columns = [acc + "_dfact_meV" for acc in accuracies]
    columns += [acc + "_ecut_deltafactor" for acc in accuracies]

    print(tabulate(data[columns], headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))
    if len(data) > 5:
        print(tabulate(data[columns].describe(), headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))

    """
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
    """

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
    """
    Check validity of pseudodojo report. Print errors to stdout
    """
    retcode = 0
    for p in options.pseudos:

        #if p.element.is_lanthanoid or p.element.is_actinoid:
        if p.element.is_actinoid:
            print("Ignoring actinoid:", os.path.relpath(p.filepath))
            continue

        rc = check_pseudo(p, check_trials=options.check_trials, verbose=options.verbose)
        retcode += rc
        if rc != 0: print(80*"=")

    if retcode != 0:
        cprint("dojo_check return code: %d [%d/%d]" % (retcode, retcode, len(options.pseudos)), "red")
    else:
        cprint("dojo_check [OK]", "green")

    return retcode


def dojo_raren(options):
    """
    Analyze results for lantanides.
    """
    from pseudo_dojo.refdata.lantanides.database import raren_database
    import pandas as pd
    retcode = 0
    with_soc = False
    icmod1, icmod3 = {}, {}
    for i, pseudo in enumerate([p for p in options.pseudos if p.element.is_lanthanoid]):
        db = raren_database(pseudo.xc)
        trial = "raren_relax" if not with_soc else "raren_relax_soc"
        if "3+" not in pseudo.basename: continue
        try:
            data = pseudo.dojo_report.get_pdframe(trial)
        except Exception:
            cprint("[%s] dojo_trial raren_relax not present!" % pseudo.basename, "red")
            continue
        ecut, a = np.array(data["ecut"])[-1], np.array(data["relaxed_a"])[-1]
        try:
            dvals = {code: db.table[code][pseudo.symbol] for code in db.table.keys()}
        except KeyError:
            cprint("Cannot find %s in pandas table!" % pseudo.symbol, "red")
            continue

        print("%s:" % pseudo.basename, "a_relax:", a, "exp: %s" % dvals["exp"])
        #for k in ("ref", "Wien2k"):
        for k in ("ref",):
            ref = dvals[k]
            if pd.isnull(ref) is None:
                cprint("Null ref for %s!" % pseudo.basename, "red")
                continue
            print("   rel_err: %.3f%% " % ( 100 * (a - ref) / ref), "wrt `%s` %.3f" % (k, ref))
        #print(dvals)
        if "icmod1" in pseudo.basename:
            icmod1[pseudo.symbol] = a
        else:
            icmod3[pseudo.symbol] = a

    table = db.table.copy()
    table["icmod1"] = [icmod1.get(s, None) for s in table.index]
    table["icmod3"] = [icmod3.get(s, None) for s in table.index]
    import matplotlib.pyplot as plt
    table[["exp", "VASP", "icmod1", "icmod3", "Wien2k"]].plot()
    plt.show()

    return retcode


@prof_main
def main():
    def str_examples():
        return """\
Usage example:
    dojodata.py plot H.psp8                ==> Plot dojo data for pseudo H.psp8
    dojodata.py compare H.psp8 H-low.psp8  ==> Plot and compare dojo data for pseudos H.psp8 and H-low.psp8
    dojodata.py nbcompare H.psp8 H-low.psp8 ==> Plot and compare dojo data in ipython notebooks.
    dojodata.py trials H.psp8 -r 1
    dojodata.py table .                    ==> Build table (find all psp8 files within current directory)
    dojodata.py figures .                  ==> Plot periodic table figures
    dojodata.py notebook H.psp8            ==> Generate ipython notebook and open it in the browser
    dojodata.py check table/*/*.psp8 -v --check-trials=gbrv_fcc,gbrv_bcc
    dojodata.py raren .
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
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('pseudos', nargs="+", help="Pseudopotential file or directory containing pseudos")
    copts_parser.add_argument('-s', "--symbols", type=parse_symbols,
        help=("List of chemical symbols to include or exclude."
              "Example --symbols=He,Li to include He and Li, --symbols=-He to exclude He"))
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='Verbose, can be supplied multiple times to increase verbosity')

    copts_parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")
    copts_parser.add_argument('--no-colors', default=False, help='Disable ASCII colors')
    copts_parser.add_argument('--seaborn', action="store_true", help="Use seaborn settings")

    # Options for pseudo selection.
    group = copts_parser.add_mutually_exclusive_group()
    group.add_argument("-r", '--rows', default="", type=parse_rows, help="Select these rows of the periodic table.")
    group.add_argument("-f", '--family', type=str, default="", help="Select this family of the periodic table.")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    plot_options_parser = argparse.ArgumentParser(add_help=False)
    plot_options_parser.add_argument("-w", "--what-plot", type=str, default="all",
                                      help="Quantity to plot e.g df for deltafactor, gbrv for GBRV tests")
    plot_options_parser.add_argument("-e", "--eos", action="store_true", help="Plot EOS curve")

    # Subparser for plot command.
    p_plot = subparsers.add_parser('plot', parents=[copts_parser, plot_options_parser],
                                   help=dojo_plot.__doc__)

    # Subparser for notebook command.
    p_notebook = subparsers.add_parser('notebook', parents=[copts_parser],
                                       help=dojo_notebook.__doc__)
    parser.add_argument('--foreground', action='store_true', default=False,
                         help="Run jupyter notebook in the foreground.")
    p_notebook.add_argument('--no-validation', action='store_true', default=False,
                             help="Don't add the validation cell.")
    p_notebook.add_argument('--hide-code', action='store_true', default=False,
                            help="Add a cell that hides the raw code.")
    p_notebook.add_argument('--no-tmp', action='store_true', default=False,
                            help="Don't use temporary file for notebook.")

    # Subparser for compare.
    p_compare = subparsers.add_parser('compare', parents=[copts_parser, plot_options_parser],
                                      help=dojo_compare.__doc__)

    # Subparser for nbcompare.
    p_nbcompare = subparsers.add_parser('nbcompare', parents=[copts_parser, plot_options_parser],
                                        help=dojo_nbcompare.__doc__)

    # Subparser for figures
    p_figures = subparsers.add_parser('figures', parents=[copts_parser], help=dojo_figures.__doc__)

    # Subparser for table command.
    p_table = subparsers.add_parser('table', parents=[copts_parser], help=dojo_table.__doc__)
    p_table.add_argument("-j", '--json', default=False, action="store_true",
                         help="Dump table in json format to file table.json")
    p_table.add_argument("-b", '--best', default=False, action="store_true",
                         help="Select best pseudos according to deltafactor")

    p_nbtable = subparsers.add_parser('nbtable', parents=[copts_parser], help=dojo_nbtable.__doc__)

    # Subparser for dist command.
    p_dist = subparsers.add_parser('dist', parents=[copts_parser], help=dojo_dist.__doc__)

    # Subparser for trials command.
    p_trials = subparsers.add_parser('trials', parents=[copts_parser], help=dojo_trials.__doc__)
    p_trials.add_argument("--savefig", type=str, default="", help="Save plot to savefig file")

    # Subparser for check command.
    def parse_trials(s):
        if s is None: return s
        #if s == "all": return DojoReport.ALL_TRIALS
        return s.split(",")

    p_check = subparsers.add_parser('check', parents=[copts_parser], help=dojo_check.__doc__)
    p_check.add_argument("--check-trials", type=parse_trials, default=None, help="List of trials to check")

    # Subparser for validate command.
    #p_validate = subparsers.add_parser('validate', parents=[copts_parser], help=dojo_validate.__doc__)

    # Subparser for raren command.
    p_raren = subparsers.add_parser('raren', parents=[copts_parser], help=dojo_raren.__doc__)

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

    if options.no_colors:
        # Disable colors
        termcolor.enable(False)

    def get_pseudos(options):
        """
        Find pseudos in paths, return :class:`DojoTable` object sorted by atomic number Z.
        Accepts filepaths or directory.
        """
        exts = ("psp8", "xml")

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
    if options.verbose: print(options.pseudos)

    if options.seaborn:
        import seaborn as sns
        sns.set(style="dark", palette="Set2")
        #sns.set(style='ticks', palette='Set2')
        #And to remove "chartjunk", do:
        #sns.despine()
        #plt.tight_layout()
        #sns.despine(offset=10, trim=True)

    # Dispatch
    return globals()["dojo_" + options.command](options)


if __name__ == "__main__":
    sys.exit(main())

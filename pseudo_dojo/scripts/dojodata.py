#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals

import sys
import os
import argparse
import pandas as pd

from warnings import warn
from collections import OrderedDict, namedtuple
from pprint import pprint
from tabulate import tabulate
from monty.os.path import find_exts
from pymatgen.io.abinit.pseudos import Pseudo
from pseudo_dojo.core.pseudos import DojoTable
from pseudo_dojo.ppcodes.oncvpsp import OncvOutputParser
from pandas import DataFrame, concat
from pymatgen.io.abinitio.netcdf import NetcdfReaderError


def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


def dojo_figures(options):
    """
    Create figures for a dojo table.
    currently for all pseudo's in the search space the one with the best df per element is chosen 
    this should probably come from a dojotable eventually
    """

    if False:
	"""
        read the data from a data file in in staed of psp files
        """
        rows = []
        with open('data') as data_file:
	    for line in data_file:
                line.rstrip('\n')
                print(line)
		data = line.split(',')
                print(data)
                data_dict = {'name': data[0],
                             'high_dfact_meV': float(data[1]),
                             'rell_high_dfact_meV': float(data[2]),
                             'high_dfactprime_meV': float(data[3])}
                if data[5] != 'nan':
                             data_dict['high_gbrv_bcc_a0_rel_err'] = float(data[5])
                             data_dict['high_gbrv_fcc_a0_rel_err'] = float(data[7])
                rows.append(data_dict)
    else:
	pseudos = options.pseudos

    	data_dojo, errors = pseudos.get_dojo_dataframe()

    	# add data that is not part of the dojo report
    	data_pseudo = DataFrame(columns=('nv', 'valence', 'rcmin', 'rcmax') )
    	for index, p in data_dojo.iterrows():
            out = p.name.replace('psp8', 'out')
	    outfile = p.symbol+'/'+out
       	    parser = OncvOutputParser(outfile)
            parser.scan()
            data_pseudo.loc[index] = [int(parser.nv), parser.valence, parser.rc_min, parser.rc_max]
  
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
                l = {k: getattr(select, k) for k in ('name', 'Z', 'high_b0_GPa', 'high_b1', 'high_v0', 'high_dfact_meV', 
                                             'high_dfactprime_meV', 'high_ecut', 'high_gbrv_bcc_a0_rel_err', 
                                             'high_gbrv_fcc_a0_rel_err', 'high_ecut', 'low_phonon', 'high_phonon',
                                             'low_ecut_hint', 'normal_ecut_hint', 'high_ecut_hint',
                                             'nv', 'valence', 'rcmin', 'rcmax')} 
            except AttributeError:
                l = {k: getattr(select, k) for k in ('name', 'Z', 'high_b0_GPa', 'high_b1', 'high_v0', 'high_dfact_meV',
                                             'high_dfactprime_meV', 'high_ecut', 
                                             'low_ecut_hint', 'normal_ecut_hint', 'high_ecut_hint',
                                             'nv', 'valence', 'rcmin', 'rcmax')}
       
            rows.append(l)

    import matplotlib.pyplot as plt
    from ptplotter.plotter import ElementDataPlotter
    import matplotlib.cm as mpl_cm
    from matplotlib.collections import PatchCollection 
    import numpy as np

    class ElementDataPlotterRangefixer(ElementDataPlotter):
        """
        modified plotter that alows to set the clim for the plot
        """
 
        def draw(self, colorbars=True, **kwargs):
            self.cbars = []
            clims = kwargs.get('clims', None)
            n = len(self.collections)
            if clims is None:
                clims = [None]*n
            elif len(clims) == 1:
                clims = [clims[0]]*n
            elif len(clims) == n:
                pass
            else:
                raise RuntimeError('incorrect number of clims provided in draw')
            for coll, cmap, label, clim  in zip(self.collections, self.cmaps, self.cbar_labels, clims):
                #print(clim)
                pc = PatchCollection(coll, cmap=cmap)
                pc.set_clim(vmin=clim[0],vmax=clim[1])
                #print(pc.get_clim())
                pc.set_array(np.array([ p.value for p in coll ]))
                self._ax.add_collection(pc)

                if colorbars:
                    options = {
                                'orientation':'horizontal',
                                'pad':0.05, 'aspect':60
                              }

                    options.update(kwargs.get('colorbar-options', {}))
                    cbar = plt.colorbar(pc, **options)
                    cbar.set_label(label)
                    self.cbars.append(cbar)
            fontdict = kwargs.get('font', {'color':'white'})
            for s in self.squares:
                if not s.label:
                    continue
                x = s.x + s.dx/2
                y = s.y + s.dy/2
                self._ax.text(x, y, s.label, ha='center',
                                             va='center',
                                             fontdict=fontdict)

            qs_labels = [k.split('[')[0] for k in self.labels]

            if self.guide_square:
                self.guide_square.set_labels(qs_labels)
                pc = PatchCollection(self.guide_square.patches, match_original=True)
                self._ax.add_collection(pc)
            self._ax.autoscale_view()

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
        try:
            return elt['high_dfact_meV']
        except KeyError:
            return float('NaN')
    def dfp(elt):
        """Delta Factor Prime"""
        try:
            return elt['high_dfactprime_meV']
        except KeyError:
            return float('NaN')
    def bcc(elt):
        """GBRV BCC [% relative error]"""
        try:
            v_bcc =  elt['high_gbrv_bcc_a0_rel_err'] if str(elt['high_gbrv_bcc_a0_rel_err']) != 'nan' else -99
        #    print(v_bcc)
            return v_bcc
        except KeyError:
            print('bcc func fail: ', elt)
            return -99 #float('NaN')

    def fcc(elt):
        """GBRV FCC [% relative error]"""
        try:
            v_fcc =  elt['high_gbrv_fcc_a0_rel_err'] if str(elt['high_gbrv_fcc_a0_rel_err']) != 'nan' else -99
        #    print(v_fcc)
            return v_fcc 
        except KeyError:
            print('fcc func fail: ', elt)
            return -99 #float('NaN')

    def low_phon_with(elt):
        """Acoustic mode low_cut """
        try:
            return elt['low_phonon'][0]
        except (KeyError, TypeError):
            #print('low_phon wiht func fail: ', elt)
            return float('NaN')

    def high_phon_with(elt):
        """AC mode [\mu eV] """
        try:
            return elt['high_phonon'][0]*1000
        except (KeyError, TypeError):
            #print('high_phon with func fail: ', elt)
            return float('NaN')
    
    def high_ecut(elt):
        """ecut high [Ha] """
        try:
            return elt['high_ecut_hint']
        except (KeyError, TypeError):
            #print('high_ecut with func fail: ', elt)
            return float('NaN')

    def low_ecut(elt):
        """ecut low [Ha] """
        try:
            return elt['low_ecut_hint']
        except (KeyError, TypeError):
            #print('low_ecut with func fail: ', elt)
            return float('NaN')

    def normal_ecut(elt):
        """ecut normal [Ha] """
        try:
            return elt['normal_ecut_hint']
        except (KeyError, TypeError):
            #print('normal_ecut with func fail: ', elt)
            return float('NaN')

    els = []
    elsgbrv = []
    elsphon = []
    rel_ers = []
    elements_data = {}    
    for el in rows:
        symbol = el['name'].split('.')[0].split('-')[0]
        try:
            rel_ers.append(max(abs(el['high_gbrv_bcc_a0_rel_err']),abs(el['high_gbrv_fcc_a0_rel_err'])))
        except (TypeError, KeyError):
            pass
        if el['high_dfact_meV'] > 0:
            elements_data[symbol] = el
            els.append(symbol)
        else:
            print('failed reading df  :', symbol, el['high_dfact_meV'])
        try:
            if el['high_gbrv_bcc_a0_rel_err'] > -100 and el['high_gbrv_fcc_a0_rel_err'] > -100:
                elsgbrv.append(symbol)
        except KeyError:
            pass
        else:
            print('failed reading gbrv: ', symbol, el['high_gbrv_bcc_a0_rel_err'], el['high_gbrv_fcc_a0_rel_err'])
            # print(el)
        try:
            if len(el['high_phonon']) > 2:
                elsphon.append(symbol)
        except (KeyError, TypeError):
            pass
      
    try:
        max_rel_err = 0.05 * int((max(rel_ers) / 0.05) + 1)
    except ValueError:
        max_rel_err = 0.20

    # plot the GBRV/DF results periodic table
    epd = ElementDataPlotterRangefixer(elements=els, data=elements_data)
    cm1 = mpl_cm.jet
    cm2 = mpl_cm.jet
    cm1.set_under('w', 1.0)
    epd.ptable(functions=[bcc,fcc,df], font={'color':color}, cmaps=[cm1,cm1,cm2], 
               clims=[[-0.6,0.6],[-0.6, 0.6], [-4,4]])
    plt.show()
    for cm2 in [mpl_cm.PiYG_r, mpl_cm.PRGn_r,mpl_cm.RdYlGn_r]:
         epd = ElementDataPlotterRangefixer(elements=els, data=elements_data)
         epd.ptable(functions=[bcc,fcc,df], font={'color':color}, cmaps=[cm1,cm1,cm2],
               clims=[[-max_rel_err,max_rel_err],[-max_rel_err, max_rel_err], [0,3]])
         plt.show()


    #plt.savefig('gbrv.eps', format='eps')

    # plot the periodic table with df and dfp
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


def dojo_plot(options):
    """Plot DOJO results for a single pseudo."""
    pseudos = options.pseudos
    if os.path.isfile('LDA'):
        xc='LDA'
    else:
        xc=None
    for pseudo in pseudos:
        if not pseudo.has_dojo_report:
            warn("pseudo %s does not contain the DOJO_REPORT section" % pseudo.filepath)
            continue

        report = pseudo.dojo_report
        #print(pseudo)
        #print(report)

        #print("trials passed: ", report.trials)
        #report.has_hints
        #report.has_exceptions
        #report.print_table()

        # ebands
        if any(k in options.what_plot for k in ("all", "ebands")):
            if report.has_trial("ebands"):
                try:
                    report.plot_ebands(title=pseudo.basename)
                except NetcdfReaderError:
                    pass

        # Deltafactor
        if report.has_trial("deltafactor") and any(k in options.what_plot for k in ("all", "df")):
            if options.eos: 
                # Plot EOS curve
                report.plot_deltafactor_eos(title=pseudo.basename)

            # Plot total energy convergence.
            fig = report.plot_deltafactor_convergence(title=pseudo.basename, xc=xc)

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


def dojo_notebook(options):
    from pseudo_dojo.pseudos import make_open_notebook
    for p in options.pseudos:
        make_open_notebook(p.filepath)


def dojo_compare(options):
    """Plot and compare DOJO results for multiple pseudos."""
    pseudos = options.pseudos
    for z in pseudos.zlist: 
        pseudos_z = pseudos[z]
        if len(pseudos_z) > 1:
            pseudos_z.dojo_compare(what=options.what_plot)
        else:
            print("Find only one pseudo for Z=%s" % z)


def dojo_trials(options):
    """Visualize the results of the different tests."""
    pseudos = options.pseudos

    # Build pandas DataFrame
    data, errors = pseudos.get_dojo_dataframe()
    #print(data)

    if errors:
        print("ERRORS:")
        pprint(errors)
    
    #import matplotlib.pyplot as plt
    #data.plot_trials(savefig=options.savefig)
    #data.plot_hist(savefig=options.savefig)
    #data.sns_plot(savefig=options.savefig)
    print(data["high_dfact_meV"])
    import matplotlib.pyplot as plt
    ax = data["high_dfact_meV"].hist()
    ax.set_xlabel("Deltafactor [meV]")
    plt.show()


def dojo_table(options):
    """Build and show a pandas table."""
    pseudos = options.pseudos

    data, errors = pseudos.get_dojo_dataframe()
    #data.tabulate()
    #return

    print(data.columns)

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

        return

    if errors:
        print("ERRORS:")
        pprint(errors)

    accuracies = ["low", "normal", "high"]
    keys = ["dfact_meV", "dfactprime_meV", "v0", "b0_GPa", "b1", "ecut_deltafactor", "ecut_hint"]
    columns = ["symbol"] + [acc + "_" + k for k in keys for acc in accuracies]
    print(columns)

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
    except:
        pass
 
    try:
        for acc in ['low', 'normal', 'high']:
            data[acc + "_abs_fcc"] = abs(data[acc + "_gbrv_fcc_a0_rel_err"])
            data[acc + "_abs_bcc"] = abs(data[acc + "_gbrv_bcc_a0_rel_err"])
    except KeyError:
        print('no GBRV data')
        pass

    wrong = data[data["high_b1"] < 0]
    if not wrong.empty:
        print("WRONG".center(80, "*") + "\n", wrong)

#    data = calc_errors(data)

    try:
        data.to_json('table.json')
    except ValueError:
        print('writing table to json failed')

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
    floatfmt=".3f"

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
    except KeyError:
        print('No GBRV data')    

    accuracies = ['low', 'normal', 'high']
    columns = [acc + "_ecut_hint" for acc in accuracies]
    print(tabulate(data[columns], headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))
    if len(data) > 5:
        print(tabulate(data[columns].describe(), headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))



    #print(data.to_string(columns=columns))

    if len(data) > 5: 
        bad = data[data["high_dfact_meV"] > data["high_dfact_meV"].mean() + data["high_dfact_meV"].std()]
        good = data[data["high_dfact_meV"] < data["high_dfact_meV"].mean()]
        print("\nPSEUDOS with high_dfact > mean plus one std (%.3f + %.3f):\n" % (data["high_dfact_meV"].mean(), data["high_dfact_meV"].std())) # ".center(80, "*"))
        print(tabulate(bad[["high_dfact_meV", "high_ecut_deltafactor"]], headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))

    #gbrv_fcc_bad = data[data["high_gbrv_fcc_a0_rerr"] > (data["high_gbrv_fcc_a0_rerr"].abs()).mean()]
    #print("\nPSEUDOS with high_dfact > mean:\n") # ".center(80, "*"))
    #print(tabulate(bad, headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))


def dojo_dist(options):
    """
    Plot the distribution of the deltafactor and of the relative error for the GBRV fcc/bcc tests.
    """
    fig = options.pseudos.plot_dfgbrv_dist()


def dojo_check(options):

    for p in options.pseudos:

        try:
            report = p.dojo_report
        except Exception as exc:
            print("Invalid dojo_report in:", p.basename)
            print("Exception: ", exc)
            continue

        #print(report)
        # Comment this to fix the md5 checksum in the pseudos
        p.check_and_fix_dojo_md5()

        #if "ppgen_hints" not in report: # and "deltafactor" not in report:
        #    print(p.basename, "old version without ppgen_hints")
        #    continue

        try:
            error = report.check()
            if error:
                print("[%s] Validation problem" % p.basename)
                print(error)
                print()

        except Exception as exc:
            print("Error: ", p.basename + str(exc))


def dojo_make_hints(options):
    from pymatgen.util.io_utils import ask_yesno, prompt
    import numpy as np

    if os.path.isfile('LDA'):
        xc='LDA'
    else:
        xc=None

    for pseudo in options.pseudos:
        try:
            report = pseudo.dojo_report

            hints = report.compute_hints()
            print("hints for %s computed from deltafactor prime: %s" % (pseudo.basename, hints))
            report.plot_deltafactor_convergence(xc=xc)           
 
            ans = ask_yesno("Do you accept the hints? [Y]")
            if ans:
                report.add_hints(hints)
                print(report.has_hints)
                pseudo.write_dojo_report(report)
            else:
                print("The dojoreport contains ecuts :\n%s" % report.ecuts)
                new_ecuts = prompt("Enter new ecuts to compute (comma-separated values or empty string to abort)")
                if len(new_ecuts) == 0:
                    print("Exit requested by user")
                    return 
                new_ecuts = np.array([float(k) for k in new_ecuts.strip().split(",")])
                print(new_ecuts)
                report.add_ecuts(new_ecuts)
                pseudo.write_dojo_report(report)

        except Exception as exc:
            print(straceback())
            print(pseudo.basename, "raised: ", str(exc))


def dojo_validate(options):
    from pymatgen.util.io_utils import ask_yesno, prompt
    from time import gmtime, strftime
    for p in options.pseudos:

        data, errors = options.pseudos.get_dojo_dataframe()

        try:
        
            #test if report is present
            try:
                report = p.dojo_report
            except Exception as exc:
                print("Invalid dojo_report in:", p.basename)
                print("Exception: ", exc)
                continue

            #test if already validated
            if 'validation' in report.keys():
                print('this pseudo was validated by %s on %s.' % (report['validation']['validated_by'], report['validation']['validated_on']))
		if not ask_yesno('Would you like to validate again? [Y]'):
                    continue

            #test if hints are present
            

            #test for ghosts
            print('\n= GHOSTS TEST ===========================================\n')
            if report.has_trial('ebands'):
                for ecut in report['ebands'].keys():
                    if "ghost_free_upto_eV" in report["ebands"][ecut].keys():
                        print('%s: Pseudo is reported to be ghost free up to %s eV' % (ecut, report["ebands"][ecut]["ghost_free_upto_eV"]))
                    else:
                        report.plot_ebands(ecut=ecut)
                        ans = float(prompt('Please enter the energy (eV) up until there is no sign of ghosts:\n'))
                        if ans > 0:
                            report["ebands"][ecut]["ghost_free_upto_eV"] = ans
                            p.write_dojo_report(report)
            else:
                print('no ebands trial present, pseudo cannot be validated')
                continue

            #test trials
            print('\n= TRIALS TEST ===========================================\n')
            try:
                error = report.check()
                if error:
                    print("[%s] Validation problem" % p.basename)
                    print(error)
                    print()

            except Exception as exc:
                print("Error: ", p.basename + str(exc))

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
                    if entry in data.keys():
                        columns.append(entry)
                        headers.append(k)
                print('ECUTS for accuracy %s:' % acc)
                print(tabulate(data[columns], headers=headers, tablefmt=tablefmt, floatfmt=floatfmt))

            #test hash

            #plot the model core charge
           
        except Exception as exc:
            print(straceback())
            print(p.basename, "raised: ", str(exc))
        
        #ask the final question
        if ask_yesno('Will you validate this pseudo? [n]'):
            name = prompt("please enter your name for later reference : ")
            report['validation'] = {'validated_by': name, 'validated_on': strftime("%Y-%m-%d %H:%M:%S", gmtime()) }
            p.write_dojo_report(report)
	

def main():
    def str_examples():
        examples = """\
    Usage example:\n
    dojodata plot H.psp8                ==> Plot dojo data for pseudo H.psp8
    dojodata trials H.psp8 -r 1
    dojodata compare H.psp8 H-low.psp8  ==> Plot and compare dojo data for pseudos H.psp8 and H-low.psp8
    dojodata table .                    ==> Build table (find all psp8 files within current directory)
    dojodata figures .                  ==> Plot periodic table figures
    dojodata notebook H.psp8            ==> Generated ipython notebook and open it in the browser
"""
        return examples

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Parent parser for commands that need to know on which subset of pseudos we have to operate.
    def parse_rows(s):
        if not s: return []
        tokens = s.split(",")
        return list(map(int, tokens)) if tokens else []

    def parse_symbols(s):
        if not s: return []
        return s.split(",")

    pseudos_selector_parser = argparse.ArgumentParser(add_help=False)
    pseudos_selector_parser.add_argument('pseudos', nargs="+", help="Pseudopotential file or directory containing pseudos")
    pseudos_selector_parser.add_argument('-s', "--symbols", type=parse_symbols, help=("List of chemical symbols to include or exclude."
        "Example --symbols=He,Li to include He and Li, --symbols=-He to exclude He"))

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
    p_plot = subparsers.add_parser('plot', parents=[pseudos_selector_parser, plot_options_parser], help="Plot DOJO_REPORT data.")

    # Subparser for notebook command.
    p_notebook = subparsers.add_parser('notebook', parents=[pseudos_selector_parser], help="Generate notebook notebook.")

    # Subparser for compare.
    p_compare = subparsers.add_parser('compare', parents=[pseudos_selector_parser, plot_options_parser], help="Compare pseudos")

    # Subparser for figures
    p_figures = subparsers.add_parser('figures', parents=[pseudos_selector_parser], help="Plot table figures")

    # Subparser for table command.
    p_table = subparsers.add_parser('table', parents=[pseudos_selector_parser], help="Build pandas table.")

    # Subparser for dist command.
    p_dist = subparsers.add_parser('dist', parents=[pseudos_selector_parser], 
                                   help="Plot distribution of deltafactor and GBRV relative errors.")

    # Subparser for trials command.
    p_trials = subparsers.add_parser('trials', parents=[pseudos_selector_parser], help="Plot DOJO trials.")
    p_trials.add_argument("--savefig", type=str, default="", help="Save plot to savefig file")

    # Subparser for check command.
    p_check = subparsers.add_parser('check', parents=[pseudos_selector_parser], help="Check pseudos")

    # Subparser for validate command.
    p_validate = subparsers.add_parser('validate', parents=[pseudos_selector_parser], help="Validate pseudos")

    # Subparser for make_hints command.
    p_make_hints = subparsers.add_parser('make_hints', parents=[pseudos_selector_parser], help="Add hints for cutoffs for pseudos")

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
        exts=("psp8",)

        paths = options.pseudos
        if len(paths) == 1 and os.path.isdir(paths[0]):
            top = os.path.abspath(paths[0])
            #print("top", top)
            paths = find_exts(top, exts, exclude_dirs="_*")
            #table = DojoTable.from_dir(paths[0])

        pseudos = []
        for p in paths:
            try:
                pseudos.append(Pseudo.from_file(p))
            except Exception as exc:
                warn("Error in %s:\n%s. This pseudo will be ignored" % (p, exc))

        table = DojoTable(pseudos)

        # Here we select a subset of pseudos according to family or rows
        if options.rows:
            table = table.select_rows(options.rows)
        elif options.family:
            table = table.select_families(options.family)

        if options.symbols:
            table = table.select_symbols(options.symbols)

        return table.sort_by_z()

    # Build DojoTable from the paths specified by the user.
    options.pseudos = get_pseudos(options)

    if options.seaborn:
        import seaborn as sns
        #sns.set(style='ticks', palette='Set2')
        sns.set(style="dark", palette="Set2")
        #And to remove "chartjunk", do:
        #sns.despine()
        #plt.tight_layout()
        #sns.despine(offset=10, trim=True)

    # Dispatch
    globals()["dojo_" + options.command](options)
    return 0


if __name__ == "__main__":
    try:
        do_prof = sys.argv[1] == "prof"
        if do_prof or do_tracemalloc: sys.argv.pop(1)
    except: 
        do_prof = False
        pass

    if do_prof:
        import pstats, cProfile
        cProfile.runctx("main()", globals(), locals(), "Profile.prof")
        s = pstats.Stats("Profile.prof")
        s.strip_dirs().sort_stats("time").print_stats()
    else:
        sys.exit(main())

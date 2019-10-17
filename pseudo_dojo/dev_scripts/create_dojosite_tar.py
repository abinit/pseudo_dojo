#!/usr/bin/env python
"""
Script to create a tar with the downloadable files and html for the pseudo-dojo.org website.

create the directory structure

element
    EX_version
    ..
    ..
..
    ..
    ..


each of these contains
    the static html version of the notebook.
    the psp8 and upf version of the pseudo
    index.html

index.html contains a list of all EX_version_type and for each entry a link to open the static html notebook and
links to download the psp8 and upf files.

the file ~/.jupyter/jupyter_notebook_config.py makes that all notebooks on save also are stored as html
http://protips.maxmasnick.com/ipython-notebooks-automatically-export-py-and-html

ipython nbconvert --to html

"""
import sys
import os
from shutil import copyfile
from pseudo_dojo.util.notebook import write_notebook_html

#PSEUDOS_TO_INCLUDE = ['ONCVPSP-PBE-PDv0.2', 'ONCVPSP-PBE-PDv0.3', 'ONCVPSP-PBEsol-PDv0.3', 'ONCVPSP-PW-PDv0.3']
PSEUDOS_TO_INCLUDE = ['ONCVPSP-PBE-PDv0.3', 'ONCVPSP-PW-PDv0.3']

#PSEUDOS_TO_INCLUDE = ['TEST', 'TEST2']

LINK_NAMES = {'djrepo': 'Dojo Report (json)', 'html': 'Dojo Report HTML', 'psp8': 'pseudo in psp8 format',
              'upf': 'pseudo in upf format'}

def main():
    website = 'website'

    usage = "Usage: " \
            "create_dojosite_tar.py path\n" \
            "creates the tar for the pseudo dojo website. To change the pseudos included change the list\n" \
            "PSEUDOS_TO_INCLUDE in create_dojosite_tar.py"
    if "--help" in sys.argv or "-h" in sys.argv:
        print(usage)
        return 1
    try:
        path = sys.argv[1]
    except:
        print(usage)
        return 1

    #  create the main directory for the website file-system

    if os.path.isdir(website):
        print('directory website already exists first move it away then rerun.')
        return
    else:
        os.makedirs(website)

    #  walk the current tree, create the directory structure and copy the .in, .psp8, and .djrepo files
    print('copying selected pseudos:\n%s' % PSEUDOS_TO_INCLUDE)

    for dirName, subdirList, fileList in os.walk(path):
        if website in dirName or 'archive' in dirName or 'DEV' in dirName:
            continue
        else:
            element = os.path.split(dirName)[-1]
            if len(element) > 2:  # in this case we did not hit an element but something else maybe better is element not in elements
                continue
            if not os.path.isdir(os.path.join(website, element)):
                os.makedirs(os.path.join(website, element))
            pseudo = os.path.split(os.path.split(dirName)[0])[1]
            if pseudo not in PSEUDOS_TO_INCLUDE:
                continue
            website_location = os.path.join(website, element, pseudo)
            if not os.path.isdir(website_location):
                os.makedirs(website_location)

        for fname in fileList:
            if fname.split('.')[-1] in ['in', 'psp8', 'djrepo', 'out']:
                copyfile(os.path.join(dirName, fname), os.path.join(website_location, fname))

    #  walk the new tree and create the .upf and notebook html files
    print('creating HTML versions of the notebooks')

    for dirName, subdirList, fileList in os.walk(website):
        for file_name in fileList:
            if 'psp8' in file_name:
                pseudopath = os.path.join(dirName, file_name)
                write_notebook_html(pseudopath)
            if '.in' in file_name:
                pass
                # todo create upf

    #  walk the new tree again and create the index.html files.
    print('creating the index HTML pages')

    HTML_string = ''  # todo add the preamble
    location = None

    for dirName, subdirList, fileList in os.walk(os.path.join(path, website)):
        if 0 < len(subdirList) < 40:  # new element todo: get rid of the magic number 40
            if location is not None:
                with open(os.path.join(location, 'index.html'), 'w') as f:
                    f.write(HTML_string)

            element = os.path.split(dirName)[-1]
            location = dirName
            HTML_string = "<H1>Pseudo's for %s</H1><BR>\n" % element
            continue

        pseudo = os.path.split(dirName)[-1]
        if len(pseudo) > 2:
            HTML_string += "<H2>%s</H2>\n<table><tr>" % pseudo
        prev_name = None
        for fname in fileList:
            if fname.split('.')[-1] in ['html', 'psp8', 'upf', 'djrepo']:
                name = fname.split('.')[0]
                if name != prev_name:
                    if prev_name is not None:
                        HTML_string += "\n"
                    HTML_string += "</tr><tr><td>%s</td>" % name
                prev_name = name
                HTML_string += "<td><A href='%s'>%s</A></td> " % (os.path.join(pseudo, fname), LINK_NAMES[fname.split('.')[-1]])
        HTML_string += "</tr></table>\n"

    with open(os.path.join(location, 'index.html'), 'w') as f:
        f.write(HTML_string)

    return

if __name__ == "__main__":
    sys.exit(main())

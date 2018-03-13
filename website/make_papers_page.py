#!/usr/bin/env python

"""
create the input for the paper page from the papers_using_pseudodojo.tex list of DOIs
"""

from operator import itemgetter
from betterbib.crossref import Crossref, pybtex_to_dict

with open('papers_using_pseudodojo.txt', 'r') as ps:
    dois = set([line.strip() for line in ps.readlines()])


def get_person_str(p):
    out = ' '.join(filter(None, [
        ' '.join(p['first'] + p['middle']),
        ' '.join(p['lineage']),
        ' '.join(p['prelast'] + p['last'],)
        ]))

    # If the name is completely capitalized, it's probably by mistake.
    if out == out.upper():
        out = out.title()
    return out


def make_ref(entry_dict):
    authors = ', '.join([get_person_str(p) for p in entry_dict['author']][:-1])
    if len(entry_dict['author']) > 2:
        authors += ','
    if len(entry_dict['author']) > 1:
        authors += ' and '
    authors += get_person_str(entry_dict['author'][-1])
    ref = '%s <BR> %s<BR>' % (entry_dict['title'], authors)
    if 'journal' in entry_dict:
        ref += '<i>%s</i>, ' % entry_dict['journal']

    if 'volume' in entry_dict:
        ref += '<b>%s</b>' % entry_dict['volume']

        if 'number' in entry_dict:
            ref += ' (%s)' % entry_dict['number']

            if 'pages' in entry_dict:
                ref += ', <i>%s</i>' % entry_dict['pages'].replace('--', '-')
    ref += ' (%s)' % entry_dict['year']
    ref += ' <a href="https://doi.org/%s">crossref</a>' % entry_dict['doi']
    return ref

html = ''

cr = Crossref()
sorted_entry_dicts = sorted([pybtex_to_dict(cr.get_by_doi(doi)) for doi in dois], key=itemgetter('year'), reverse=True)

for entry_dict in sorted_entry_dicts:
    html += '<blockquote>' + make_ref(entry_dict) + '</blockquote>'

with open('papers.html', 'w') as page:
    page.write(html)

#!/usr/bin/env python

"""
create the input for the paper page from the papers_using_pseudodojo.tex list of DOIs
"""
import json
import pandas
from operator import itemgetter
from betterbib.crossref import Crossref, pybtex_to_dict

with open('papers_using_pseudodojo.txt', 'r') as ps:
    dois = [doi for doi in set([line.strip() for line in ps.readlines()]) if '#' not in doi]


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
    #print("entry_dict", entry_dict)
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
            else:
                ref += ', <i>%s</i>' % entry_dict['doi'].split('.')[-1]

    ref += ' (%s)' % entry_dict['year']
    ref += ' <a href="https://doi.org/%s">crossref</a>' % entry_dict['doi']
    return ref


def make_tweet(doi):
    cr = Crossref()
    entry = pybtex_to_dict(cr.get_by_doi(doi))

    tweet = f'New #compchem work using #PseudoDojo pseudopotentials: {entry["title"]} by {entry["author"][0]["last"][0]} et al. https://doi.org/{doi}'

    return tweet


cr = Crossref()
sorted_entry_dicts = sorted([pybtex_to_dict(cr.get_by_doi(doi)) for doi in dois], key=itemgetter('year'), reverse=True)
df = pandas.DataFrame(sorted_entry_dicts)
ax = df['year'].groupby(df["year"]).count().plot(kind="bar")
ax.set_ylabel('publications per year')
fig = ax.get_figure()
fig.set_size_inches(8, 4)
fig.tight_layout()
fig.savefig('year.png')
ax.get_figure().savefig('year.png', transparent=True)

print(json.dumps(sorted_entry_dicts, indent=2))

t1 = '<table><tr><td width="100" align="center" valign="top"><div data-badge-popover="right" ' \
     'class="altmetric-embed" data-badge-type="donut" ' \
     'data-doi="'
t2 = '" /></td><td>'
t3 = '</td></tr></table><BR />'

html = ''.join([t1 + entry_dict['doi'] + t2 + make_ref(entry_dict) + t3 for entry_dict in sorted_entry_dicts])

with open('papers.html', 'w') as page:
    page.write(html)

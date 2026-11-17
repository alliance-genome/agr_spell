"""Dataset metadata from NCBI: GEO series (title, summary, PubMed IDs, experiment type) and PubMed.

Uses the public E-utilities (esearch and esummary), in batches and at most about three
requests a second, as NCBI asks of clients without an API key.
"""
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
BATCH = 150
PAUSE = 0.4
TIMEOUT = 120
GSE = re.compile(r'GSE\d+')


class NCBIError(Exception):
    pass


def _post(tool, **params):
    """One E-utilities call (POST, so long ID lists fit); JSON, or XML text for efetch. Retries twice."""
    params.setdefault('retmode', 'json')
    params.update(tool='spell-next')
    data = urllib.parse.urlencode(params).encode()
    for attempt in range(3):
        try:
            with urllib.request.urlopen(EUTILS + tool, data=data, timeout=TIMEOUT) as resp:
                out = json.load(resp) if params['retmode'] == 'json' else resp.read().decode('utf-8', 'replace')
            time.sleep(PAUSE)
            return out
        except (OSError, ValueError) as err:
            if attempt == 2:
                raise NCBIError('%s failed: %s' % (tool, err))
            time.sleep(3)


def geo_id_of(filename):
    """The GEO series in a file name (GSE123456.final.pcl -> GSE123456), or None."""
    m = GSE.search(filename)
    return m.group(0) if m else None


def assay_of(gdstype):
    """microarray / RNA-seq from GEO's experiment type, or None when it says neither."""
    gdstype = (gdstype or '').lower()
    if 'by high throughput sequencing' in gdstype:
        return 'RNA-seq'
    if 'by array' in gdstype:
        return 'microarray'
    return None


def fetch_series(gse_ids, post=_post):
    """{GSE: {title, summary, pubmed_ids, gdstype}} for the series GEO knows."""
    out = {}
    ids = sorted(set(gse_ids))
    for i in range(0, len(ids), BATCH):
        chunk = ids[i:i + BATCH]
        term = '(%s) AND gse[ETYP]' % ' OR '.join('%s[ACCN]' % g for g in chunk)
        uids = post('esearch.fcgi', db='gds', term=term, retmax=len(chunk) * 2)['esearchresult']['idlist']
        if not uids:
            continue
        result = post('esummary.fcgi', db='gds', id=','.join(uids))['result']
        for uid in result.get('uids', []):
            s = result[uid]
            if s.get('accession') in chunk:
                out[s['accession']] = {'title': s.get('title') or '', 'summary': s.get('summary') or '',
                                       'pubmed_ids': [str(p) for p in s.get('pubmedids') or []],
                                       'gdstype': s.get('gdstype') or ''}
    return out


def fetch_pubmed(pmids, post=_post):
    """{PMID: {title, authors, first_author, journal, year}} in SGD's style ("Gasch AP", "Cell")."""
    out = {}
    ids = sorted(set(pmids))
    for i in range(0, len(ids), BATCH):
        result = post('esummary.fcgi', db='pubmed', id=','.join(ids[i:i + BATCH]))['result']
        for uid in result.get('uids', []):
            s = result[uid]
            authors = [a['name'] for a in s.get('authors', []) if a.get('authtype', 'Author') == 'Author']
            year = re.match(r'\d{4}', s.get('pubdate') or '')
            out[uid] = {'title': s.get('title') or '', 'authors': ', '.join(authors),
                        'first_author': authors[0] if authors else s.get('sortfirstauthor') or '',
                        'journal': s.get('fulljournalname') or s.get('source') or '',
                        'year': year.group(0) if year else ''}
    return out


def fetch_abstracts(pmids, post=_post):
    """{PMID: abstract text} (labelled sections joined, e.g. "BACKGROUND: ... RESULTS: ...")."""
    out = {}
    ids = sorted(set(pmids))
    for i in range(0, len(ids), BATCH):
        root = ET.fromstring(post('efetch.fcgi', db='pubmed', id=','.join(ids[i:i + BATCH]), retmode='xml'))
        for art in root.iter('PubmedArticle'):
            pmid = art.findtext('.//MedlineCitation/PMID')
            parts = []
            for t in art.findall('.//Abstract/AbstractText'):
                text = ' '.join(''.join(t.itertext()).split())
                parts.append('%s: %s' % (t.get('Label'), text) if t.get('Label') else text)
            if pmid and parts:
                out[pmid] = ' '.join(parts)
    return out


def metadata_for(gse_ids, post=_post):
    """{GSE: dataset table fields} for each series: publication details when GEO links a PubMed record
    (title, authors, journal, year and abstract from PubMed), else the GEO series title and summary;
    assay from GEO's experiment type."""
    series = fetch_series(gse_ids, post)
    pmids = [s['pubmed_ids'][0] for s in series.values() if s['pubmed_ids']]
    papers = fetch_pubmed(pmids, post)
    abstracts = fetch_abstracts(pmids, post) if pmids else {}
    out = {}
    for gse, s in series.items():
        paper = papers.get(s['pubmed_ids'][0]) if s['pubmed_ids'] else None
        out[gse] = {
            'geo_id': gse,
            'pubmed_id': s['pubmed_ids'][0] if paper else '',
            'title': (paper or {}).get('title') or s['title'],
            'description': (abstracts.get(s['pubmed_ids'][0]) if paper else None) or s['summary'],
            'first_author': (paper or {}).get('first_author', ''),
            'authors': (paper or {}).get('authors', ''),
            'journal': (paper or {}).get('journal', ''),
            'year': (paper or {}).get('year', ''),
            'assay': assay_of(s['gdstype']) or '',
            'gdstype': s['gdstype'],
        }
    return out

#!/usr/bin/env python3
"""Turn a WormBase SPELL release folder into a spell-next raw directory.

WormBase publishes each release as a folder (e.g. TablesForSPELL_9.3.0) generated from ACeDB:
gene_list.txt, dataset_list.txt, noempty/ (PCL files), dataset_table_enriched.txt (the old
SPELL datasets table), systematic_to_common.txt, alias_to_systematic.txt, gene_association.wb
and gene_ontology.obo.

    python3 tools/wb_import.py SRC_DIR RAW_DIR [--zeros-as-missing]

writes RAW_DIR/gene_list.txt, dataset_list.txt and expression/ (copies of the listed PCL files),
and RAW_DIR/meta/genes.tsv (gene_table format), datasets.tsv (dataset_table format),
gene_ontology.obo and gene_association.wb. mod.toml is not touched (start from configs/wb.toml).

Dataset tags in the enriched table look like "Method: RNAseq|Species: Brugia malayi|Topic: diet":
Method becomes the assay and Topic and "Tissue Specific" entries become tags. Each gene's species
comes from WBGeneSpe.txt (WBGene ID, species) when the release has it, else from the species of the
datasets measuring it. HTML in titles is removed.
WormBase writes zero expression as log2(1e-10) = -33.2. Such values dominate correlations, so
--zeros-as-missing writes them as empty cells (not measured) in the copied PCL files.
Uses only the Python standard library, so it runs on the server without the spell-next image.
"""
import csv
import html
import os
import re
import shutil
import sys
from collections import Counter, defaultdict

ZERO = -33.2  # WormBase writes log2(1e-10) = -33.219 for zero expression


def clean(text):
    """Plain text from an HTML fragment (tags dropped, entities decoded, whitespace collapsed)."""
    return ' '.join(html.unescape(re.sub(r'<[^>]+>', ' ', text or '')).split())


def links_of(field):
    """Links from the enriched title field: 'Title <br> GEO Record: <a ...>GSE1</a> Platform: <a>GPL2</a>
    <br> Lab Atlas: <a>...</a> <br>Download ...: <a>file.csv</a>' -> 'GEO Record GSE1 <url> | ...'."""
    out = []
    for section in re.split(r'<br\s*/?>', field or '', flags=re.I)[1:]:
        prefix, pos = '', 0
        for m in re.finditer(r'<a[^>]+href="([^"]+)"[^>]*>(.*?)</a>', section, flags=re.I | re.S):
            before = clean(section[pos:m.start()]).strip(' :,')
            if before:                          # a new label such as "GEO Record" or "Platform"
                prefix = 'Download' if before.lower().startswith('download') else before
            pos = m.end()
            text = clean(m.group(2))
            label = prefix or text if text.startswith(('http://', 'https://')) else ('%s %s' % (prefix, text)).strip()
            out.append('%s <%s>' % (label.replace('|', '/').replace('<', '').replace('>', ''), m.group(1).strip()))
    return ' | '.join(out)


def parse_tags(field):
    """'Method: RNAseq|Species: X|Topic: a|Topic: b' -> (method, species, [topics])."""
    method, species, topics = '', '', []
    for part in (field or '').split('|'):
        key, _, value = part.partition(':')
        key, value = key.strip().lower(), value.strip()
        if not value:
            continue
        if key == 'method':
            method = value
        elif key == 'species':
            species = value
        elif key in ('topic', 'tissue specific'):
            topics.append(value)
    return method, species, topics


def read_table(path):
    """Rows of the enriched datasets table: id, pmid, filename, GDS, GPL, channels, title, description,
    conditions, genes, first author, authors, publication title, journal, year, condition links, tags."""
    with open(path, encoding='utf-8', errors='replace') as fh:
        return [line.rstrip('\n').split('\t') for line in fh if line.strip()]


def main(src, raw, zeros_as_missing=False):
    os.makedirs(os.path.join(raw, 'expression'), exist_ok=True)
    os.makedirs(os.path.join(raw, 'meta'), exist_ok=True)
    for name in ('gene_list.txt', 'dataset_list.txt'):
        shutil.copy(os.path.join(src, name), os.path.join(raw, name))
    for name in ('gene_ontology.obo', 'gene_association.wb'):
        shutil.copy(os.path.join(src, name), os.path.join(raw, 'meta', name))

    listed = []
    with open(os.path.join(src, 'dataset_list.txt')) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                listed.append(parts[1])
    table = {r[2]: r for r in read_table(os.path.join(src, 'dataset_table_enriched.txt')) if len(r) >= 17}

    # Datasets: the listed ones are searched; the other table rows are shown for reference only
    # (display_only), e.g. single-cell datasets with no expression matrix.
    cols = ['dataset_id', 'filename', 'geo_id', 'pubmed_id', 'title', 'description', 'first_author', 'authors',
            'journal', 'year', 'tags', 'assay', 'value_type', 'species', 'display_only', 'links']
    species_of_dataset, assays, missing = {}, Counter(), []
    listed_set = set(listed)
    display = [n for n in table if n not in listed_set]
    with open(os.path.join(raw, 'meta', 'datasets.tsv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(cols)
        for name in listed + display:
            r = table.get(name)
            if r is None:
                missing.append(name)
                w.writerow(['', name] + [''] * (len(cols) - 2))
                continue
            method, species, topics = parse_tags(r[16])
            species_of_dataset[name] = species
            if name in listed_set:
                assays[method] += 1
            geo = re.search(r'GEO Record:.*?\b(GSE\d+)\b', r[6])
            title = clean(r[12]) or clean(r[6].split('<br>')[0])
            w.writerow([r[0] if r[0].isdigit() else '', name, geo.group(1) if geo else '',
                        r[1] if r[1].isdigit() else '', title, clean(r[7]), clean(r[10]), clean(r[11]),
                        clean(r[13]), r[14] if r[14].isdigit() else '',
                        '|'.join(dict.fromkeys(t.replace('|', '/') for t in topics)), method, '', species,
                        '' if name in listed_set else 'yes', links_of(r[6])])

    # Expression files, and each gene's species from the datasets measuring it.
    gene_species = defaultdict(Counter)
    zero_files, zero_values, values = 0, 0, 0
    for name in listed:
        species = species_of_dataset.get(name, '')
        zeros = 0
        with open(os.path.join(src, 'noempty', name), encoding='utf-8', errors='replace') as fh, \
                open(os.path.join(raw, 'expression', name), 'w', encoding='utf-8') as out:
            for i, line in enumerate(fh):
                f = line.rstrip('\n').split('\t')
                if i == 0 or f[0].upper() == 'EWEIGHT' or not f[0]:
                    out.write(line)
                    continue
                if species:
                    gene_species[f[0]][species] += 1
                for k in range(3, len(f)):
                    try:
                        x = float(f[k])
                    except ValueError:
                        continue
                    values += 1
                    if x <= ZERO:
                        zeros += 1
                        if zeros_as_missing:
                            f[k] = ''
                out.write('\t'.join(f) + '\n' if zeros_as_missing else line)
        zero_values += zeros
        zero_files += zeros > 0

    # Genes: WBGene ID is both the systematic name (PCL rows, gene_list.txt) and the GAF ID.
    common = {}
    with open(os.path.join(src, 'systematic_to_common.txt'), encoding='utf-8', errors='replace') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 3:
                common[f[1]] = f[2]
    aliases = defaultdict(dict)
    with open(os.path.join(src, 'alias_to_systematic.txt'), encoding='utf-8', errors='replace') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            # Skip cross-reference strings such as CAEEL/WormBase=WBGene.../UniProtKB=...
            if len(f) >= 3 and f[1] and f[1] != common.get(f[2]) and '=' not in f[1]:
                aliases[f[2]][f[1]] = None
    listed_species = {}
    spe_path = os.path.join(src, 'WBGeneSpe.txt')
    if os.path.exists(spe_path):
        with open(spe_path, encoding='utf-8', errors='replace') as fh:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) >= 2 and f[0].startswith('WBGene') and f[1]:
                    listed_species[f[0]] = f[1]
    genes, no_species, named, disagree = 0, 0, 0, 0
    with open(os.path.join(raw, 'gene_list.txt')) as gl, \
            open(os.path.join(raw, 'meta', 'genes.tsv'), 'w', newline='', encoding='utf-8') as out:
        w = csv.writer(out, delimiter='\t', lineterminator='\n')
        w.writerow(['systematic', 'mod_id', 'standard', 'aliases', 'species', 'description'])
        for line in gl:
            f = line.rstrip('\n').split('\t')
            if len(f) < 2:
                continue
            wb = f[1]
            measured = gene_species[wb].most_common(1)[0][0] if gene_species.get(wb) else ''
            sp = listed_species.get(wb) or measured
            disagree += bool(listed_species.get(wb) and measured and measured != listed_species[wb])
            no_species += not sp
            named += wb in common
            w.writerow([wb, wb, common.get(wb, ''), '|'.join(a.replace('|', '/') for a in aliases.get(wb, {})),
                        sp, ''])
            genes += 1

    print('datasets: %d listed (searched), %d display-only (shown without expression data); %d listed '
          'without a table row' % (len(listed), len(display), len(missing)))
    print('assays: ' + ', '.join('%s %d' % kv for kv in assays.most_common()))
    print('genes: %d (%d with a common name, %d with aliases); species from %s; %d without a species%s'
          % (genes, named, len(aliases), 'WBGeneSpe.txt (%d genes) and datasets' % len(listed_species)
             if listed_species else 'the datasets measuring them', no_species,
             '; %d measured in datasets of another species' % disagree if disagree else ''))
    print('zero placeholders (<= %.1f, i.e. log2 of ~0): %d of %d values (%.2f%%), in %d of %d files%s'
          % (ZERO, zero_values, values, 100.0 * zero_values / max(values, 1), zero_files, len(listed),
             '; written as missing' if zeros_as_missing else ''))


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    if len(args) != 2 or set(sys.argv[1:]) - set(args) - {'--zeros-as-missing'}:
        sys.exit(__doc__)
    main(args[0], args[1], zeros_as_missing='--zeros-as-missing' in sys.argv)

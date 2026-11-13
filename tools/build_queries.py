"""Build a validation query set from real SPELL searches in the Rails log extract.

Usage: build_queries.py RAILS_TSV_GZ GENE_LIST OUT_TSV
Output columns: query_id, kind (single|multi), comma-separated systematic gene names.
"""
import gzip
import random
import re
import sys
from collections import Counter

SEP = re.compile(r'[\s,;]+')


def main(rails_tsv, gene_list, out_path):
    names = set()
    with open(gene_list) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) == 2:
                names.add(parts[1].upper())

    singles, multis = Counter(), Counter()
    with gzip.open(rails_tsv, 'rt', errors='replace') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 8 or f[2] != 'SearchController#show_results' or not f[7]:
                continue
            genes = tuple(dict.fromkeys(g.upper() for g in SEP.split(f[7].strip()) if g))
            if not genes or not all(g in names for g in genes):
                continue
            if len(genes) == 1:
                singles[genes] += 1
            elif len(genes) <= 20:
                multis[genes] += 1

    rng = random.Random(42)
    top = [g for g, _ in singles.most_common(20)]
    rest = [g for g in singles if g not in top]
    chosen_single = top + rng.sample(rest, 30)
    multi_pool = sorted(multis)
    chosen_multi = rng.sample(multi_pool, min(30, len(multi_pool)))

    with open(out_path, 'w') as out:
        for i, genes in enumerate(chosen_single, 1):
            out.write('s%02d\tsingle\t%s\n' % (i, ','.join(genes)))
        for i, genes in enumerate(chosen_multi, 1):
            out.write('m%02d\tmulti\t%s\n' % (i, ','.join(genes)))
    print('distinct resolvable single-gene queries: %d, multi-gene: %d' % (len(singles), len(multis)))
    print('wrote %d single + %d multi queries to %s' % (len(chosen_single), len(chosen_multi), out_path))


if __name__ == '__main__':
    main(*sys.argv[1:4])

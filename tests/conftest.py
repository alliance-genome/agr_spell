"""Shared test data: a small two-species MOD written in the raw-data layout."""
import gzip

import numpy as np

SPECIES = {'A': 'Alpha worm', 'B': 'Beta worm'}
GENES_PER_SPECIES = 12
N_CONDS = 10

CONFIG = '''
[mod]
id = "TM"
name = "Test MOD"
organism = "test"
examples = ["a-1", "a-1 a-2"]

[links]
gene = "https://example.org/gene/{mod_id}"

[genes]
file = "meta/genes.tsv"
format = "gene_table"

[datasets]
metadata = "meta/datasets.tsv"
format = "dataset_table"

[go]
obo = "meta/go.obo"
gaf = ["meta/a.gaf.gz", "meta/b.gaf"]
'''

OBO = '''format-version: 1.2
data-version: releases/test

[Term]
id: GO:0000001
name: everything
namespace: biological_process

[Term]
id: GO:0000002
name: alpha module
namespace: biological_process
is_a: GO:0000001 ! everything

[Term]
id: GO:0000003
name: beta module
namespace: biological_process
is_a: GO:0000001 ! everything
'''


def gaf_line(mod_id, term):
    return '\t'.join(['TM', mod_id, mod_id, '', term, 'PMID:1', 'IDA', '', 'P', '', '', 'gene', 'taxon:1',
                      '20260101', 'TM', '', '']) + '\n'


def make_toy_raw(raw, n_datasets_per_species=2, seed=7):
    """Write a two-species MOD (gene table, dataset table, PCLs, GO files, mod.toml) into raw.

    Each species has an RNA-seq dataset (set0, tagged stress) and a microarray one (set1, tagged growth).
    """
    (raw / 'expression').mkdir(parents=True, exist_ok=True)
    (raw / 'meta').mkdir(parents=True, exist_ok=True)
    (raw / 'mod.toml').write_text(CONFIG)
    (raw / 'meta' / 'go.obo').write_text(OBO)

    gene_list, gene_rows = [], ['systematic\tmod_id\tstandard\taliases\tspecies\tdescription']
    gid = 1000
    for sp, species in SPECIES.items():
        for i in range(1, GENES_PER_SPECIES + 1):
            gene_list.append('%d\t%s%d' % (gid, sp, i))
            gene_rows.append('%s%d\tMOD%s%d\t%s-%d\t%s%d-alias\t%s\tgene %d of %s'
                             % (sp, i, sp, i, sp.lower(), i, sp.lower(), i, species, i, species))
            gid += 1
    (raw / 'gene_list.txt').write_text('\n'.join(gene_list) + '\n')
    (raw / 'meta' / 'genes.tsv').write_text('\n'.join(gene_rows) + '\n')

    # Two datasets per species; genes 1-4 of each species rise and fall together.
    rng = np.random.default_rng(seed)
    ds_rows, ds_list = ['filename\tgeo_id\ttitle\tfirst_author\tyear\ttags\tassay\tvalue_type'], []
    ds_id = 1
    for sp in SPECIES:
        for k in range(n_datasets_per_species):
            name = '%s_set%d.pcl' % (sp, k)
            base = rng.normal(size=N_CONDS)
            lines = ['ID\tNAME\tGWEIGHT\t' + '\t'.join('c%d' % c for c in range(N_CONDS)),
                     'EWEIGHT\t\t\t' + '\t'.join(['1'] * N_CONDS)]
            for i in range(1, GENES_PER_SPECIES + 1):
                vals = (base + 0.1 * rng.normal(size=N_CONDS)) if i <= 4 else rng.normal(size=N_CONDS)
                lines.append('%s%d\t%s%d\t1\t' % (sp, i, sp, i) + '\t'.join('%.4f' % (8 + v) for v in vals))
            (raw / 'expression' / name).write_text('\n'.join(lines) + '\n')
            (raw / 'expression' / (name + '.sinfo')).write_text('All\t0\t1\n')
            ds_list.append('%d\t%s' % (ds_id, name))
            ds_rows.append('%s\tGSE%d\t%s study %d\tAuthor %s\t2020\t%s\t%s' % (
                name, 9000 + ds_id, SPECIES[sp], k, sp, 'stress' if k == 0 else 'growth',
                'RNA-seq\tlog2 TPM' if k == 0 else 'microarray\t'))
            ds_id += 1
    (raw / 'dataset_list.txt').write_text('\n'.join(ds_list) + '\n')
    (raw / 'meta' / 'datasets.tsv').write_text('\n'.join(ds_rows) + '\n')

    # Species A's co-expressed genes share "alpha module"; species B's share "beta module".
    with gzip.open(raw / 'meta' / 'a.gaf.gz', 'wt') as fh:
        fh.write('!date-generated: 20260101\n')
        fh.writelines(gaf_line('MODA%d' % i, 'GO:0000002' if i <= 4 else 'GO:0000001')
                      for i in range(1, GENES_PER_SPECIES + 1))
    (raw / 'meta' / 'b.gaf').write_text('!date-generated: 20260102\n' + ''.join(
        gaf_line('MODB%d' % i, 'GO:0000003' if i <= 4 else 'GO:0000001') for i in range(1, GENES_PER_SPECIES + 1)))

    return raw

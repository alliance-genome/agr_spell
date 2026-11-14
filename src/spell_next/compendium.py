"""Load SPELL's expression compendium (gene list, dataset list, PCL and .sinfo files).

Parsing follows Spell++ (ExpressionDataset.cpp, utilities.cpp) so that search results
can be compared with the production engine. With ``faithful=True`` two Spell++ quirks
are reproduced:

* runs of tabs are collapsed when splitting a line, so an empty cell shifts the values
  after it one column to the left;
* reading stops at the first blank line.
"""
import math
import os
import re
from dataclasses import dataclass, field

import numpy as np

_FLOAT_PREFIX = re.compile(r'\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)')


def tokenize(line, faithful=True):
    """Split a tab-delimited line the way utilities::tokenize does (or plainly)."""
    if faithful:
        return [tok for tok in line.split('\t') if tok != '']
    return line.split('\t')


def parse_float(token):
    """Mimic ``istringstream >> float``: parse the longest numeric prefix, else NaN."""
    match = _FLOAT_PREFIX.match(token)
    return float(match.group(1)) if match else math.nan


def read_gene_list(path):
    """Return (ids, names) from SPELL's gene_list.txt (``<id>\\t<systematic name>``)."""
    ids, names = [], []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2 and parts[0].strip():
                ids.append(int(parts[0]))
                names.append(parts[1].strip())
    return ids, names


def read_dataset_list(path):
    """Return [(dataset_id, filename)] from SPELL's dataset_list.txt."""
    out = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2 and parts[0].strip():
                out.append((int(parts[0]), parts[1].strip()))
    return out


def read_sinfo(path):
    """Return (z_mean, z_std) from the ``All`` line of a .sinfo stats file."""
    with open(path) as fh:
        parts = fh.readline().rstrip('\n').split('\t')
    return parse_float(parts[1]), parse_float(parts[2])


@dataclass
class Dataset:
    """One expression dataset: a gene x condition matrix plus its z statistics."""

    id: int
    filename: str
    conditions: list
    gene_index: np.ndarray          # row -> index into Compendium.gene_ids
    values: np.ndarray              # float64, rows x conditions, NaN = missing
    z_mean: float = math.nan
    z_std: float = math.nan
    row_of: dict = field(default_factory=dict, repr=False)  # gene index -> row

    @property
    def num_conds(self):
        return len(self.conditions)


def read_pcl(path, name_to_index, faithful=True):
    """Parse a PCL file into (conditions, gene_index array, values matrix)."""
    with open(path, errors='replace') as fh:
        lines = fh.read().split('\n')
    header = tokenize(lines[0].rstrip('\r'), faithful)
    offset = 3 if len(header) > 2 and header[2].upper() == 'GWEIGHT' else 2
    conditions = header[offset:]
    num_conds = len(conditions)

    pos = 1
    if pos < len(lines) and tokenize(lines[pos], faithful)[:1] and \
            tokenize(lines[pos], faithful)[0].upper() == 'EWEIGHT':
        pos += 1

    rows = {}
    for line in lines[pos:]:
        if line == '' and faithful:
            break
        tokens = tokenize(line, faithful)
        if not tokens:
            continue
        gene = None
        for cand in (tokens[0], tokens[0].upper(),
                     tokens[1] if len(tokens) > 1 else '',
                     tokens[1].upper() if len(tokens) > 1 else ''):
            if cand in name_to_index:
                gene = name_to_index[cand]
                break
        if gene is None:
            continue
        vals = [parse_float(t) for t in tokens[offset:]]
        if len(vals) > num_conds:
            raise ValueError('too many conditions for %s in %s' % (tokens[0], path))
        vals += [math.nan] * (num_conds - len(vals))
        rows[gene] = vals           # a later duplicate row replaces an earlier one

    gene_index = np.fromiter(rows.keys(), dtype=np.int32, count=len(rows))
    values = np.array(list(rows.values()), dtype=np.float64).reshape(len(rows), num_conds)
    return conditions, gene_index, values


class Compendium:
    """All genes and datasets SPELL searches over."""

    def __init__(self, gene_ids, gene_names, datasets):
        self.gene_ids = list(gene_ids)
        self.gene_names = list(gene_names)
        self.index_of_id = {gid: i for i, gid in enumerate(self.gene_ids)}
        self.index_of_name = {name.upper(): i for i, name in enumerate(self.gene_names)}
        self.datasets = datasets
        for ds in datasets:
            ds.row_of = {int(g): r for r, g in enumerate(ds.gene_index)}

    @classmethod
    def from_spell_files(cls, raw_dir, faithful=True, dedupe=False):
        """Load from a copy of SPELL's gene_list.txt, dataset_list.txt and expression/.

        ``dedupe=True`` keeps only the first dataset ID listed for each file.
        """
        gene_ids, gene_names = read_gene_list(os.path.join(raw_dir, 'gene_list.txt'))
        name_to_index = {name: i for i, name in enumerate(gene_names)}
        expr_dir = os.path.join(raw_dir, 'expression')
        parsed = {}
        datasets = []
        for ds_id, filename in read_dataset_list(os.path.join(raw_dir, 'dataset_list.txt')):
            if dedupe and filename in parsed:
                continue
            if filename not in parsed:
                path = os.path.join(expr_dir, filename)
                conds, gene_index, values = read_pcl(path, name_to_index, faithful)
                # .sinfo holds Spell++'s z statistics; new datasets need none (the build recomputes them).
                sinfo = path + '.sinfo'
                z_mean, z_std = read_sinfo(sinfo) if os.path.exists(sinfo) else (math.nan, math.nan)
                parsed[filename] = (conds, gene_index, values, z_mean, z_std)
            conds, gene_index, values, z_mean, z_std = parsed[filename]
            datasets.append(Dataset(ds_id, filename, conds, gene_index, values, z_mean, z_std))
        return cls(gene_ids, gene_names, datasets)

    def save(self, path):
        """Cache the parsed compendium as a single .npz file."""
        arrays = {
            'gene_ids': np.array(self.gene_ids, dtype=np.int64),
            'gene_names': np.array(self.gene_names),
            'ds_ids': np.array([d.id for d in self.datasets], dtype=np.int64),
            'ds_files': np.array([d.filename for d in self.datasets]),
            'ds_z': np.array([[d.z_mean, d.z_std] for d in self.datasets]),
        }
        for k, d in enumerate(self.datasets):
            arrays['g%d' % k] = d.gene_index
            arrays['v%d' % k] = d.values
            arrays['c%d' % k] = np.array(d.conditions, dtype=object).astype(str)
        np.savez(path, **arrays)

    @classmethod
    def load(cls, path):
        """Load a compendium cached by :meth:`save`."""
        with np.load(path) as npz:
            datasets = []
            for k, (ds_id, fname) in enumerate(zip(npz['ds_ids'], npz['ds_files'])):
                z_mean, z_std = npz['ds_z'][k]
                datasets.append(Dataset(int(ds_id), str(fname), list(npz['c%d' % k]),
                                        npz['g%d' % k], npz['v%d' % k],
                                        float(z_mean), float(z_std)))
            return cls(npz['gene_ids'].tolist(), npz['gene_names'].tolist(), datasets)

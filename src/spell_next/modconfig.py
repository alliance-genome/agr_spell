"""Per-MOD configuration (one config file per SPELL instance: SGD, WormBase, ...).

The build reads ``mod.toml`` from the raw data directory; paths in it are relative to
that directory. Without a config file the SGD defaults below apply, which match the
layout copied from spell.yeastgenome.org. The resolved config is stored in the build
(``build_info``), so the web app needs only the build output.

Example for a multi-species MOD (see configs/ for complete files)::

    [mod]
    id = "WB"
    name = "WormBase"
    organism = "nematodes"
    contact = "help@wormbase.org"
    examples = ["unc-54", "lin-12 glp-1"]

    [links]
    gene = "https://wormbase.org/species/all/gene/{mod_id}"

    [genes]
    file = "meta/genes.tsv"
    format = "gene_table"

    [go]
    gaf = ["meta/c_elegans.gaf.gz"]
"""
import copy
import os
import tomllib

GENE_FORMATS = ('sgd_features', 'gene_table')
DATASET_FORMATS = ('spell_dump', 'dataset_table')
LINK_FIELDS = ('mod_id', 'systematic', 'standard')

DEFAULTS = {
    'mod': {
        'id': 'SGD',
        'name': 'Saccharomyces Genome Database',
        'site_title': 'SPELL next',
        'organism': 'yeast',
        'contact': 'sgd-helpdesk@lists.stanford.edu',
        'examples': ['ACT1', 'GAL1', 'HSP104', 'RPL3 RPL4A RPS2 RPS3', 'CDC28 CLN1 CLN2 CLB5 CLB6'],
        'listing_example': 'GAL1 GAL7 GAL10',
    },
    'links': {
        'gene': 'https://www.yeastgenome.org/locus/{mod_id}',
    },
    'genes': {
        'file': 'meta/SGD_features.tab',
        'format': 'sgd_features',
        'default_species': 'Saccharomyces cerevisiae',
    },
    'datasets': {
        'metadata': None,               # None: the newest meta/SQLdump* file, if any
        'format': 'spell_dump',
        'assays': ['microarray', 'RNA-seq'],   # allowed values of the dataset table's assay column
    },
    'go': {
        'obo': 'meta/go-basic.obo',
        'gaf': ['meta/gene_association.sgd.gaf.gz'],
    },
    'updates': {},                      # {raw-relative path: download URL}, fetched by `spell update-refs`
}


class ConfigError(ValueError):
    pass


def _merge(base, override):
    out = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _merge(out[key], value)
        else:
            out[key] = value
    return out


def load(raw_dir, filename='mod.toml'):
    """Return the resolved config for raw_dir (defaults merged with mod.toml if present)."""
    path = os.path.join(raw_dir, filename)
    user = {}
    if os.path.exists(path):
        with open(path, 'rb') as fh:
            user = tomllib.load(fh)
    unknown = set(user) - set(DEFAULTS)
    if unknown:
        raise ConfigError('unknown section(s) in %s: %s' % (path, ', '.join(sorted(unknown))))
    cfg = _merge(DEFAULTS, user)
    if isinstance(cfg['go']['gaf'], str):
        cfg['go']['gaf'] = [cfg['go']['gaf']]
    if cfg['datasets']['metadata'] is None and cfg['datasets']['format'] == 'spell_dump':
        meta = os.path.join(raw_dir, 'meta')
        dumps = sorted(f for f in os.listdir(meta) if f.startswith('SQLdump')) if os.path.isdir(meta) else []
        cfg['datasets']['metadata'] = os.path.join('meta', dumps[-1]) if dumps else None
    validate(cfg)
    return cfg


def validate(cfg):
    if cfg['genes']['format'] not in GENE_FORMATS:
        raise ConfigError('genes.format must be one of %s' % ', '.join(GENE_FORMATS))
    if cfg['datasets']['format'] not in DATASET_FORMATS:
        raise ConfigError('datasets.format must be one of %s' % ', '.join(DATASET_FORMATS))
    template = cfg['links']['gene']
    try:
        template.format(**{f: 'x' for f in LINK_FIELDS})
    except (KeyError, IndexError) as err:
        raise ConfigError('links.gene may only use {%s}: %s' % ('}, {'.join(LINK_FIELDS), err))
    if not isinstance(cfg['mod']['examples'], list):
        raise ConfigError('mod.examples must be a list of queries')
    assays = cfg['datasets']['assays']
    if not isinstance(assays, list) or not all(isinstance(a, str) and a for a in assays):
        raise ConfigError('datasets.assays must be a list of assay names')
    for rel, url in cfg['updates'].items():
        if os.path.isabs(rel) or '..' in rel.replace('\\', '/').split('/'):
            raise ConfigError('updates: %r must be a path inside the raw directory' % rel)
        if not isinstance(url, str) or not url.startswith(('https://', 'http://', 'file://')):
            raise ConfigError('updates: %r must map to an http(s) URL' % rel)


def path(raw_dir, relative):
    """Resolve a config path relative to the raw data directory (None stays None)."""
    return None if relative is None else os.path.join(raw_dir, relative)


def public(cfg):
    """The part of the config the web app needs (no file paths)."""
    return {'mod': cfg['mod'], 'links': cfg['links'], 'assays': cfg['datasets']['assays']}

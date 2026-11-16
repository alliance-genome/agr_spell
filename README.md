# spell-next

Prototype replacement for SGD's SPELL co-expression search, written in Python.
`spell_next.search` ports the Spell++ algorithm (ZSCORE2 distance, no orthology,
as run on spell.yeastgenome.org).

## Setup

```bash
python3.13 -m venv .venv
.venv/bin/pip install -e '.[dev]'
.venv/bin/pytest
```

## Data

Copy the live compendium from the SPELL server into `data/raw/` (not committed):

- `gene_list.txt`, `dataset_list.txt` from `/var/www/html/spell_web/config/`
- each listed `*.pcl` and `*.pcl.sinfo` from `/var/www/html/data/expression/`

Then cache it:

```python
from spell_next.compendium import Compendium
Compendium.from_spell_files('data/raw', faithful=True).save('data/compendium_faithful.npz')
```

## Curator workflow

Each instance lives in one directory (`ROOT`, or set `SPELL_ROOT`):

```
ROOT/raw/       inputs curators edit: mod.toml, gene_list.txt, dataset_list.txt, expression/, meta/
ROOT/builds/    one directory per build (YYYYMMDD-HHMMSS); the last 3 are kept
ROOT/current    symlink to the build being served
ROOT/logs/      monthly search logs (kept across builds)
```

| Command | What it does |
| --- | --- |
| `spell init ROOT --config configs/sgd.toml` | Create the layout and copy a config template |
| `spell validate ROOT` | Check inputs without building; errors block a rebuild, warnings do not; writes `validation_report.tsv` |
| `spell rebuild ROOT` | Validate, build into `builds/<time>`, smoke-test, compare with the live build, switch `current`, prune |
| `spell rebuild ROOT --no-switch` | Build and compare, but keep serving the current build |
| `spell rollback ROOT [--to BUILD]` | Serve the previous (or a chosen) build again |
| `spell status ROOT` | Live build and build history |
| `spell serve ROOT --port 8000` | Run the web app on the live build |
| `spell migrate-metadata ROOT` | One-off: turn the old SPELL metadata export into `meta/datasets.tsv` |

`rebuild` and `rollback` run `--restart-cmd` (or `$SPELL_RESTART_CMD`) afterwards; a running
server keeps serving the build it started with until it is restarted.

**Adding a dataset:** put the PCL file in `raw/expression/`, add `<new ID>\t<file>` to
`raw/dataset_list.txt` (`spell validate` prints the next free ID), add a row to
`raw/meta/datasets.tsv`, then `spell validate ROOT` and `spell rebuild ROOT`. No `.sinfo`
file is needed. To remove a dataset, delete its line from `dataset_list.txt` and rebuild.

Validation on the current SGD inputs: 0 errors, 265 warnings, including 99 datasets with
fewer than 5 conditions (never searched, in old SPELL either), 56 where under half the
rows match known genes, and 24 without metadata.

## Deployment

`Dockerfile` builds one image for all MODs (AG Grid included, checksum-pinned).
`deploy/docker-compose.yml` runs one container per MOD and `deploy/spellctl` wraps the
curator commands; see [deploy/README.md](deploy/README.md).

## One instance per MOD

Each MOD (SGD, WormBase, ...) runs its own instance from the same code, configured by
`mod.toml` in its raw data directory (paths relative to that directory). Without it,
SGD defaults apply. Templates: `configs/sgd.toml`, `configs/wb.example.toml`.

| Section | What it sets |
| --- | --- |
| `[mod]` | ID, name, organism, contact, example queries shown on the home page |
| `[links]` | Gene page URL template: `{mod_id}`, `{systematic}`, `{standard}` |
| `[genes]` | Gene names file: `sgd_features` (SGD_features.tab) or `gene_table` |
| `[datasets]` | Dataset metadata: `spell_dump` (old SPELL export) or `dataset_table` |
| `[go]` | GO ontology and one or more GAF files (matched on GAF column 2 = `mod_id`) |

`gene_table` is a tab-separated file with a header row: `systematic` (required; the
name used in the PCL files), `mod_id`, `standard`, `aliases` (`|`-separated),
`species`, `description`, `qualifier`. `dataset_table` likewise: `filename` (required),
`geo_id`, `pubmed_id`, `title`, `description`, `first_author`, `authors`, `journal`,
`year`, `tags` (`|`-separated), `assay` (one of `[datasets] assays`, by default
`microarray` or `RNA-seq`) and `value_type` (what the values are, e.g. `log2 TMM CPM`,
`log2 TPM`, `vst`, `log2 ratio`). The site filters searches, listings and the dataset list
by assay; `spell validate` warns about missing or unknown assays, RNA-seq datasets without a
`value_type`, and RNA-seq values that are not log-scale (raw counts or TPM/FPKM).

Several species can share one instance (e.g. WormBase's nematodes). Each dataset's
species is the species of most of its genes; a gene only occurs in its own species'
datasets, so results stay within the query's species. GO enrichment tests a gene list
against the annotated genes of its own species. Species columns and filters appear in
the interface only when an instance has more than one species. The resolved config is
stored in the build, so the web app needs only the build directory.

## Build and run the prototype

For SGD, `data/raw/meta/` also needs `SGD_features.tab` (from downloads.yeastgenome.org)
and SPELL's dataset metadata export (`SQLdump*_fordownloadfiles.txt` from the server);
copy `configs/sgd.toml` to `data/raw/mod.toml`.

```bash
.venv/bin/python -m spell_next.build data/raw data/build      # about 4 minutes on 16 cores
SPELL_NEXT_DATA=data/build .venv/bin/uvicorn spell_next.api:app --port 8765
```

Open http://127.0.0.1:8765. The build uses the corrected calculations (64-bit, no tab
collapsing, duplicate datasets dropped) and precomputes the top 100 genes and datasets
for every gene. Multi-gene searches run live (at most 4 at a time).

| Endpoint | Purpose |
| --- | --- |
| `GET /api/genes/lookup?q=ACT1 GAL1` | Resolve systematic, standard or alias names |
| `GET /api/genes/{id}/coexpression` | Precomputed single-gene result |
| `POST /api/search` `{"genes": [...]}` | Live multi-gene search |
| `GET /api/datasets/{id}` | Dataset metadata and conditions |
| `GET /api/datasets/{id}/expression?genes=...` | Gene-centred values for heatmaps |
| `POST /api/go-enrichment` `{"gene_ids": [...], "aspect": "P"}` | GO term enrichment |
| `GET /api/datasets` | All datasets with publication details, tags and QC flags |
| `GET /api/tags` | Dataset tags with counts |
| `GET /api/expression-listing?genes=...&tags=...` | Datasets measuring the genes (listing mode, no search) |
| `GET /api/datasets/{id}/download` | Expression matrix as searched, tab-separated |
| `GET /api/info` | Counts and data versions |
| `GET /api/config` | MOD branding, gene link template, examples and species |
| `GET /embed?gene=ACT1&genes=10&sets=5` | Compact result page for iframes |

`POST /api/search` also takes `"tags": [...]`: the search then runs only over datasets
carrying any of those tags (old SPELL searched every dataset and only hid the others).

Pages: search results (paged AG Grid tables with CSV download, GO enrichment, heatmaps of
all listed datasets with scale and colour options, click a row for the gene's full
profile), `#/datasets` (tag filter, text filter, search within tags), `#/expression/...`
(listing mode), dataset pages (QC, download), `#/help` and `#/about`.

Searches are logged monthly to `SPELL_NEXT_LOG_DIR` (default `<data>/logs`;
empty disables) as `date, time, kind, genes, tags`, with no IP addresses.

GO enrichment replaces SPELL's GO::TermFinder service. It needs `go-basic.obo`
(purl.obolibrary.org/obo/go/go-basic.obo) and `gene_association.sgd.gaf.gz`
(downloads.yeastgenome.org/curation/literature/) in `data/raw/meta/`. Annotations are
propagated over is_a and part_of; each term with at least two list genes gets a
hypergeometric test against all genes annotated in that aspect, with Bonferroni and
Benjamini-Hochberg corrections. The results page tests the top 50 hits (query genes
excluded) and shows terms with FDR < 5%. Ontology and annotation versions are stored
in the `build_info` table and shown under the results.
`spell update-refs` refreshes these files (and `SGD_features.tab`) from the URLs under
`[updates]` in `mod.toml`; on the server this runs weekly (see `deploy/README.md`).

## Data QC

The build runs `spell_next.qc` on every dataset (skip with `--no-qc`):

- unlogged intensities or counts (4 datasets) get log2(x+1); unlogged ratios (7) get log2;
- datasets of unclear scale (2) are left as is and flagged;
- condition shifts above 1 log unit (100 datasets) are flagged, not corrected;
- each dataset's z_mean/z_std are recomputed from the data actually searched.

Results go to the `dataset_qc` table, `qc_report.tsv` (for curators) and dataset pages.

Centring each condition on its median is available (`--centre`) but off, because it
made results worse. Genome-wide benchmark (`tools/benchmark_go.py`, `tools/benchmark_ci.py`):
share of a gene's top 20 hits that share a specific GO process term (5-100 genes),
over 4,565 annotated genes, with paired bootstrap 95% CIs against corrected maths without QC:

| Version | Score | Difference |
| --- | --- | --- |
| Corrected maths, no QC | 0.1118 | — |
| + log transforms (default) | 0.1117 | −0.0000 (−0.0004 to +0.0003) |
| + condition centring | 0.1080 | −0.0037 (−0.0048 to −0.0027) |
| Old SPELL arithmetic (float32) | 0.1140 | +0.0023 (+0.0015 to +0.0030) |
| Old SPELL parsing, float64 | — | +0.0002 (−0.0002 to +0.0005) |
| Current: + centred genes, capped Fisher input | 0.1143 | +0.0025 (+0.0019 to +0.0032) |
| Random hits | 0.0103 | |

The first float64 port still had two numerical problems, fixed in the current version:
near-constant genes gave correlations of ±1 through rounding (genes are now centred
before correlating), and |r| was capped at 0.999999 only when it reached 1 exactly, so
values a rounding step below 1 gave Fisher z ~19 and let one degenerate gene pair (zero
in all conditions but one) outweigh hundreds of datasets (|r| is now always capped).
With both fixed, the float64 version scores above old SPELL; the float32 version's
earlier edge came from the same kind of instability, not better biology.

## Build speed

`spell_next.precompute` computes all single-gene results dataset by dataset (each
dataset's correlations once, in row blocks) instead of one search per gene; it matches
the per-gene search exactly (`tests/test_precompute.py`; `build --per-gene` runs the slow
reference). QC computes each dataset's pair statistics once unless normalization changed
it. BLAS runs single-threaded in each worker process (multi-threaded BLAS in 3 workers
oversubscribed the container's CPU quota). A full SGD build with 3 workers takes ~3.5
minutes on an Apple M-series Mac and ~19 minutes on the m6i.xlarge dev instance (was 78).

## Validation against production

1. `tools/build_queries.py` picks real single- and multi-gene searches from the Rails log.
2. `tools/capture_baseline.sh` sends them one at a time to the live Spell++ engine
   (validated IDs only; it stops at the first error).
3. `tools/compare.py data/compendium_faithful.npz data/validation/query_ids.tsv
   data/validation/baseline.tsv --precision=float32` compares rankings.

With `--precision=float32` and `faithful=True`, all 80 validation queries match
production exactly (gene scores, dataset order and weights, up to ties).

## Known Spell++ defects (reproduced only in faithful mode)

| Defect | Effect | Correct mode |
| --- | --- | --- |
| One-pass Pearson formula in 32-bit floats | Some correlations are badly wrong (0.857 computed as 1.000); the #1 hit changes in 13 of 80 test queries | `precision='float64'` |
| Runs of tabs collapsed when reading PCL files | Values after an empty cell shift left; 25 of 802 files affected, 2 entirely | `faithful=False` |
| `stddev()` returns the variance | Only orders datasets for single-gene queries | kept as is |
| Two PCL files listed twice in `dataset_list.txt` | GSE113118 and GSE137261 counted twice | fix the list |

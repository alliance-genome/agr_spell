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

## Build and run the prototype

`data/raw/meta/` also needs `SGD_features.tab` (from downloads.yeastgenome.org) and
SPELL's dataset metadata export (`SQLdump*_fordownloadfiles.txt` from the server).

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
| Random hits | 0.0103 | |

The float32 edge comes from the unstable correlation formula, which inflates
correlations among abundant, low-variance genes; those are heavily co-annotated in GO,
so the benchmark rewards the error. Correct arithmetic is kept.

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

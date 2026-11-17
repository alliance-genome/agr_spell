# SPELL curator guide

How to add, change and remove datasets in SPELL (spell-next), and how to put the changes
on the site. Written for SGD; WormBase works the same way with `wb` in place of `sgd`
and `/srv/spell/wb` in place of `/srv/spell/sgd`.

Site: https://spell.dev.yeastgenome.org/

## Contents

1. [The short version](#1-the-short-version)
2. [Logging in and where things are](#2-logging-in-and-where-things-are)
3. [Adding datasets](#3-adding-datasets)
4. [Editing dataset metadata and tags](#4-editing-dataset-metadata-and-tags)
5. [Retiring and restoring datasets](#5-retiring-and-restoring-datasets)
6. [Validate, rebuild, check, roll back](#6-validate-rebuild-check-roll-back)
7. [What the validate messages mean](#7-what-the-validate-messages-mean)
8. [Gene names and GO: the weekly update](#8-gene-names-and-go-the-weekly-update)
9. [Troubleshooting](#9-troubleshooting)

## 1. The short version

Every change follows the same four steps:

1. **Edit the input files** under `/srv/spell/sgd/raw/` (expression files, `dataset_list.txt`,
   `meta/datasets.tsv`).
2. **Validate**: `./spellctl validate sgd`. This only checks; it changes nothing. Fix any errors.
3. **Rebuild**: `./spellctl rebuild sgd`. About 20 minutes, in the background. The site keeps
   showing the old version until the new one is ready, then switches by itself.
4. **Check** the site. If something is wrong, `./spellctl rollback sgd` puts the previous
   version back immediately.

Run `spellctl` from `~/spell-next/deploy`. Several changes can go into one rebuild.

## 2. Logging in and where things are

Log in with the instance's SSH key (ask Shuai for the key file and address):

```bash
ssh -i agr-sgd-dev-ec2.pem ubuntu@<instance>
cd ~/spell-next/deploy
```

To copy files from your computer to the server, use `scp` with the same key, e.g.
`scp -i agr-sgd-dev-ec2.pem GSE123456.final.pcl ubuntu@<instance>:/srv/spell/sgd/raw/expression/`.
Editors on the server: `emacs` or `nano`.

Everything SPELL serves is built from the files in `/srv/spell/sgd/raw/`:

| Path | What it is | Who edits it |
| --- | --- | --- |
| `raw/expression/` | One expression file (PCL) per dataset | Curators |
| `raw/dataset_list.txt` | Dataset ID and file name for every dataset SPELL builds | Curators |
| `raw/meta/datasets.tsv` | Title, authors, tags, assay, ... for each dataset | Curators |
| `raw/meta/keywords.csv` | SGD keywords: categories and definitions for tags | Curators |
| `raw/gene_list.txt` | SPELL's gene IDs and systematic names | Rarely (new genes) |
| `raw/meta/SGD_features.tab` | Gene names, aliases, descriptions | Automatic, weekly |
| `raw/meta/go-basic.obo`, `gene_association.sgd.gaf.gz` | GO ontology and annotations | Automatic, weekly |
| `raw/mod.toml` | Instance settings (site name, file locations, download URLs) | Rarely |
| `raw/retired/` | Datasets taken out of SPELL, kept for restoring | Curators |

Other locations in `/srv/spell/sgd/`:

| Path | What it is |
| --- | --- |
| `builds/` | The last 3 builds; `current` points to the live one |
| `validation_report.tsv` | Per-dataset results of the last validate (opens in Excel) |
| `current/qc_report.tsv` | Data quality results for each dataset in the live build |
| `rebuild.log`, `update.log` | Output of the last rebuild and of the weekly updates |
| `logs/searches_YYYY_MM.tsv` | What was searched each month (no IP addresses) |

**Keep a backup before editing a file by hand**, e.g.
`cp -p meta/datasets.tsv meta/datasets.tsv.bak-$(date +%Y%m%d)`.

## 3. Adding datasets

### 3.1 The expression file (PCL)

Tab-separated text, one row per gene, one column per sample:

```
YORF       NAME    GWEIGHT  heat 0 min  heat 5 min  heat 10 min  ...
EWEIGHT                     1           1           1            ...
YAL001C    TFC3    1        0.12        -0.40       0.88         ...
```

- **Rows** are keyed by systematic name (e.g. `YAL001C`) in the first column. Rows whose name
  is not in `raw/gene_list.txt` are ignored (validate reports the share that matched).
- **The `EWEIGHT` row is optional.**
- **Columns** are samples or conditions. Give each a clear, unique name; these are shown on the
  dataset page and under the heatmaps.
- **Missing values**: leave the cell empty.
- **At least 5 columns.** SPELL does not search datasets with fewer than 5 conditions, so we no
  longer load them.
- **One value type per file.** Don't mix, for example, counts and TPM in one file.
- **A `.sinfo` file is not needed**; the build computes those statistics itself.

### 3.2 Which values: log scale

SPELL compares genes within each dataset, so every file must be on a log scale and normalised
the same way across all of its columns. Units may differ between datasets.

| Data | Supply | `value_type` |
| --- | --- | --- |
| Microarray, two-colour | log2 ratios | `log2 ratio` |
| Microarray, one-colour | log2 intensities | `log2 intensity` |
| RNA-seq counts (raw, RSEM, DESeq2-normalised, CPM) | log2(TMM-normalised CPM + 1) | `log2 TMM CPM` |
| RNA-seq, TPM / FPKM / RPKM | log2(TPM + 1); convert FPKM/RPKM to TPM first | `log2 TPM` |
| RNA-seq already variance-stabilised | as supplied | `vst` or `rlog` |
| RNA-seq log ratios to a reference (e.g. time 0) | as supplied | `log2 ratio` |

If a file is not on a log scale, the build applies log2 (or log2(x + 1)) itself and validate
says so. For RNA-seq that is not good enough (raw counts are not normalised for library size),
so validate warns: fix the file rather than relying on it.

### 3.3 Step by step: `spellctl add-datasets`

1. **Copy the expression files** to `/srv/spell/sgd/raw/expression/`. Put the GEO series in the
   name, e.g. `GSE123456.final.pcl`.
2. **Optionally, prepare a metadata table** for them and put it under `/srv/spell/sgd/`, e.g.
   `/srv/spell/sgd/raw/new_datasets.csv`. Comma- or tab-separated, with a header row: `filename`
   plus any `datasets.tsv` columns (section 4). Typically `filename, geo_id, pubmed_id, title,
   description, first_author, authors, journal, year, tags, assay, value_type`.
3. **Preview**, then add:
   ```bash
   ./spellctl add-datasets sgd --metadata /srv/spell/sgd/raw/new_datasets.csv --dry-run
   ./spellctl add-datasets sgd --metadata /srv/spell/sgd/raw/new_datasets.csv
   ```
   For each file it assigns the next free dataset ID (never reusing a retired one), appends
   a line to `dataset_list.txt` and a row to `datasets.tsv`, and keeps backups of both
   (`*.bak-<date>`).
   - **Your table's values are used as given.** Anything it leaves blank (title, PubMed ID,
     authors, journal, year, abstract, assay) is filled from GEO and PubMed, using the GSE ID in
     the table or the file name. Datasets your table completes are not looked up at all.
   - **Without `--metadata`** it adds every file in `expression/` that isn't listed yet, with
     metadata from GEO and PubMed. `--tags`, `--assay` and `--value-type` set those columns for
     all of them; file names can be given to add only some.
   - **Duplicates are checked first.** A file whose data is already in SPELL (identical values,
     under any name) is refused, listed as `DUPLICATE`; use `--allow-duplicates` only if that is
     really intended. Files with the same condition names or the same GEO series as an existing
     dataset are listed as `note:` and added (a series split into several datasets is normal).
     The first run reads every dataset (under a minute); after that only new files are read.
   - **Filled-in metadata is a draft.** GEO may link several papers (the first is used), and
     SGD often gives datasets split from one series their own titles. Check the new rows.
4. **Check the new rows** in `raw/meta/datasets.tsv`: tags, assay, value_type, titles.
5. **Validate**, then **rebuild** (section 6). Add many datasets and rebuild once.

To check the datasets already in SPELL for duplicates, run `./spellctl find-duplicates sgd`. It
lists groups with the same data, the same condition names, or the same GEO series.

Adding by hand works too: put the ID and file name in `dataset_list.txt` (`./spellctl validate sgd`
prints the next free ID) and a row in `datasets.tsv`.

## 4. Editing dataset metadata and tags

`raw/meta/datasets.tsv` holds one row per dataset, tab-separated, with a header row:

| Column | Content |
| --- | --- |
| `filename` | Must match the name in `dataset_list.txt` (required) |
| `geo_id` | e.g. `GSE123456` |
| `pubmed_id` | Number only |
| `title` | Shown everywhere. A row without a title counts as "no metadata" |
| `description` | Summary shown on the dataset page |
| `first_author` | e.g. `Gasch AP` |
| `authors` | Full author list |
| `journal`, `year` | |
| `tags` | SGD keywords, separated by `\|`, e.g. `stress\|heat shock` |
| `assay` | `microarray` or `RNA-seq` |
| `value_type` | What the values are (section 3.2) |

- **Tags describe the biology** and must be SGD keywords, spelled exactly as in the keyword list
  (e.g. `QTL`, not `QTLs`). Don't use tags for the assay; that is what the `assay` column is for.
- **Assay** is required: datasets without one are left out when users filter by assay.
- **Editing in Excel**: save as "Tab delimited text" and keep the header row. Avoid line breaks
  inside cells.
- **If the metadata comes from the SGD database**, replace `datasets.tsv` with the new export
  (same columns), and make any one-off fixes in the database too, or the next export will undo
  them.

Metadata changes need a rebuild like any other change.

### Keywords and their categories

`raw/meta/keywords.csv` is the list of SGD keywords that tags may use, exported from the SGD
database. Columns: `category`, `display_name`, `format_name`, `obj_url`, `description`.

- **The tag picker groups tags by `category`**, in the order categories first appear in the file.
  Categories only group keywords for display: datasets keep their keywords, and searches still
  filter by keyword.
- **`description`** is shown when hovering over a tag; **`format_name`** builds the link to the
  keyword's SGD page from each dataset page.
- **To move a keyword to another category** (or rename or add a category), edit the `category`
  column, validate and rebuild. The file can be edited in Excel; save it as CSV.
- **Validate warns about any dataset tag that is not in this file** (a typo, or a keyword not
  yet exported), and about keywords with no category; those appear under "Other" on the site.

## 5. Retiring and restoring datasets

To take a dataset out of SPELL without losing it:

1. Move its file (and `.sinfo`, if any) from `raw/expression/` to `raw/retired/expression/`.
2. Move its line from `raw/dataset_list.txt` to `raw/retired/dataset_list.txt`.
3. Move its row from `raw/meta/datasets.tsv` to `raw/retired/datasets.tsv`.
4. Validate and rebuild.

To restore it, move the three pieces back and rebuild. In November 2026 the 99 datasets with
fewer than 5 conditions were retired this way.

## 6. Validate, rebuild, check, roll back

Run these in `~/spell-next/deploy`:

| Task | Command |
| --- | --- |
| Check inputs (changes nothing) | `./spellctl validate sgd` |
| Build and go live | `./spellctl rebuild sgd` |
| Follow a running rebuild (Ctrl-C stops following, not the rebuild) | `./spellctl rebuild-log sgd` |
| Live build, build history, rebuild in progress | `./spellctl status sgd` |
| Undo the last rebuild | `./spellctl rollback sgd` |
| Go back to a specific build | `./spellctl rollback sgd 20261116-194648` |

About a rebuild:

- **It validates first** and stops on errors, leaving the live site unchanged.
- **It runs in the background**, so you can log out. Only one rebuild runs at a time.
- **It takes about 20 minutes** for SGD: about 7 for data checks, 10 to precompute every gene.
- **Before switching, it tests the new build** and compares it with the live one ("0 datasets
  added, 99 removed; 0 genes have more than 2 of their top-10 hits changed"). A big change in
  top hits after only a metadata edit means something went wrong.
- **The last 3 builds are kept**, so a rollback is instant.

After a rebuild, check:

1. The new or changed dataset on the Datasets page and its own page (title, tags, assay,
   conditions).
2. Its "Data quality" section on the dataset page, or `current/qc_report.tsv`.
3. A search for a gene you know well.

## 7. What the validate messages mean

**Errors** stop a rebuild:

| Message | Fix |
| --- | --- |
| `listed in dataset_list.txt but not in expression/` | Copy the file, or fix the name |
| `dataset ID N used more than once` | Give one of them a new ID |
| `cannot be read: ...` | Usually a row with more values than there are column names |
| `no rows match a gene in gene_list.txt` | Wrong first column (not systematic names), or wrong organism |
| `no numeric values` | Empty or non-numeric file |

**Warnings** don't stop a rebuild. Ones to act on:

| Message | Meaning |
| --- | --- |
| `N conditions: datasets with fewer than 5 are not searched` | Retire it or don't load it |
| `RNA-seq values look like raw read counts` | Convert to log2(TMM CPM + 1) |
| `RNA-seq values are not log-scale` | Convert to log2(TPM + 1) or log2(TMM CPM + 1) |
| `N RNA-seq datasets have no value_type` | Fill in `value_type` |
| `N datasets have no assay` | Fill in `assay` |
| `assay 'X' is not one of microarray, RNA-seq` | Typo in `assay` |
| `N datasets have no metadata title` | Add the row or title in `datasets.tsv` |
| `N metadata rows name files not in dataset_list.txt` | Leftover rows, or a file name mismatch |
| `tag 'X' is not in the keyword list` | Typo in `tags`, or add the keyword to `keywords.csv` |
| `N keywords have no category` | Fill in `category` in `keywords.csv` |
| `only N% of rows match known genes` | Check the gene names in the first column |
| `duplicate condition names`, `empty condition name(s)` | Rename the columns |

Warnings that are usually fine:

| Message | Meaning |
| --- | --- |
| `will be log2 transform of unlogged ...` (info) | The file isn't logged; the build logs it. Fine for microarray |
| `conditions shifted by up to N log units` | Columns have different medians; usually fine, over ~3 worth a look |
| `unclear scale ...: left as is` | The build couldn't tell whether the values are logged; check the file |
| `listed twice (IDs A and B)` | Two IDs for one file; only the first is built |

The full list for every dataset is in `/srv/spell/sgd/validation_report.tsv`.

## 8. Gene names and GO: the weekly update

Gene names, aliases and descriptions (`SGD_features.tab`) and GO files are downloaded from SGD and
the GO site every **Sunday at 10:00 UTC**. If anything changed, SPELL rebuilds by itself. If a
download looks wrong (empty, an error page, much smaller than before), it is not used and the
current file stays. If a rebuild is already running, that week is skipped.

- To update now (e.g. after an urgent gene rename): `./spellctl update sgd`
- To see what happened: `./spellctl status sgd`, or `/srv/spell/sgd/update.log`

**New genes**: expression rows are matched against `raw/gene_list.txt`, not against
`SGD_features.tab`. A gene new to SGD needs a line in `gene_list.txt` (an unused numeric ID, e.g. the largest plus one,
a tab, its systematic name) before its expression values are used.

## 9. Troubleshooting

| Problem | What to do |
| --- | --- |
| Rebuild failed | `./spellctl rebuild-log sgd` shows why. The old build is still live. Fix and rebuild |
| New build looks wrong | `./spellctl rollback sgd`, then investigate |
| Site down or not responding | `./spellctl status sgd`; then `./spellctl restart sgd` (see `~/Restart_app.README`) |
| "A rebuild of sgd is already running" | Wait for it: `./spellctl rebuild-log sgd` |
| Changes don't show on the site | Check the rebuild finished (`./spellctl status sgd`); reload with Cmd+Shift+R |

Never run `git pull` or `./spellctl image` while a rebuild is running.

More detail for server administration: `~/spell-next/deploy/README.md` and `~/Restart_app.README`.
For anything else, contact Shuai.

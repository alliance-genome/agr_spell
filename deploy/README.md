# Running spell-next with Docker

One image serves every MOD; each MOD runs in its own container with its own data
directory on the host. `spellctl` wraps the `spell` curator commands so they run in a
one-off container and restart the app afterwards.

| MOD | Container | Host data directory | Port (localhost only) |
| --- | --- | --- | --- |
| SGD | `spell-sgd` | `/srv/spell/sgd` | 8001 |
| WormBase | `spell-wb` (profile `wb`) | `/srv/spell/wb` | 8002 |

Each data directory holds `raw/` (inputs curators edit), `builds/`, `current` (the
build being served) and `logs/` (monthly search logs). Override the location with
`SPELL_DATA_ROOT`.

## Opening the site

SGD dev: **https://spell.dev.yeastgenome.org/** (public). How it is wired:

| Piece | Setting |
| --- | --- |
| DNS | CNAME `spell.dev.yeastgenome.org` → `app-qa-lb` (Route 53 zone yeastgenome.org, sgd-prod account) |
| HTTPS | `app-qa-lb` (agr-sgd-dev), certificate `*.dev.yeastgenome.org`; HTTP redirects to HTTPS |
| Routing | listener rule 75: host `spell.dev.yeastgenome.org` → target group `spell-dev` (instance port 8001, health check `/api/info`) |
| Instance | `deploy/.env` sets `SPELL_BIND=0.0.0.0`; security group `spell-dev-from-alb` allows 8001 only from the load balancer |

Without a load balancer, ports stay on localhost. With the instance's SSH key:

```bash
ssh -i agr-sgd-dev-ec2.pem -L 8001:localhost:8001 ubuntu@<instance>
# then open http://localhost:8001 (SGD) or forward 8002 for WormBase
```

## Curator tasks (on the instance, in ~/spell-next/deploy)

| Task | Command |
| --- | --- |
| Check inputs after editing `raw/` | `./spellctl validate sgd` |
| Build and go live (validates first; keeps the old build) | `./spellctl rebuild sgd` |
| Follow a rebuild (or see how the last one ended) | `./spellctl rebuild-log sgd` |
| Undo the last rebuild | `./spellctl rollback sgd` |
| Go back to a specific build | `./spellctl rollback sgd 20261113-160000` |
| See builds and container state | `./spellctl status sgd` |
| Follow the app log | `./spellctl logs sgd` |

### Adding a dataset

`./spellctl add-datasets sgd [--metadata TABLE]` does steps 2 and 3 below, filling metadata from
a curator's table and from GEO/PubMed; see `docs/CURATOR_GUIDE.md`.

All of a MOD's inputs live under `/srv/spell/<mod>/raw/` (SGD shown; use `wb` for WormBase).

1. **Copy the expression file** to `raw/expression/`:

   ```bash
   scp GSE123456.final.pcl spell-dev:/srv/spell/sgd/raw/expression/
   ```

   PCL format: tab-separated, header `YORF  NAME  GWEIGHT  <condition 1>  <condition 2> ...`,
   an optional `EWEIGHT` row, then one row per gene keyed by systematic name (e.g.
   `YFL039C`). Leave missing values empty. log2 data is expected; unlogged intensities or
   ratios are detected, converted and flagged by QC. Datasets with fewer than 5 conditions
   are listed but not searched. A `.sinfo` file is optional (the build computes it).

2. **Register it** in `raw/dataset_list.txt`: append `<ID><TAB><file name>`, e.g.
   `3051	GSE123456.final.pcl`. `./spellctl validate sgd` prints the next free ID. IDs must be
   unique but need not be in order.

3. **Add its metadata** to `raw/meta/datasets.tsv`, one tab-separated row:

   ```
   filename  geo_id  pubmed_id  title  description  first_author  authors  journal  year  tags  assay  value_type
   ```

   `filename` must match `dataset_list.txt`. `tags` are `|`-separated (e.g.
   `stress|heat shock`) and describe the biology; reuse existing tag names so the datasets
   page filter groups them. `assay` is `microarray` or `RNA-seq` (not a tag; validate warns
   if it is missing). For RNA-seq, supply log2(TMM CPM + 1) or log2(TPM + 1) values and say
   which in `value_type` (e.g. `log2 TMM CPM`, `log2 TPM`, `vst`, `log2 ratio`); validate
   warns about raw counts or other unlogged values. A row without a `title` counts as no
   metadata (validate warns).

4. **Validate** (changes nothing): `./spellctl validate sgd`. It checks the PCL format,
   unknown genes, duplicate IDs and missing metadata; fix errors before continuing.

5. **Rebuild**: `./spellctl rebuild sgd`, then `./spellctl rebuild-log sgd` or
   `./spellctl status sgd` to follow it (see below).

6. **Check** the dataset on the site's Datasets page and its QC flags in
   `/srv/spell/sgd/current/qc_report.tsv`. If anything is wrong, `./spellctl rollback sgd`
   returns to the previous build instantly (the last 3 builds are kept).

Add several datasets at once and rebuild once; the rebuild time hardly changes. To remove a
dataset, delete its line from `dataset_list.txt` (and optionally its metadata row) and
rebuild. Files must be readable by uid 1000 (the `ubuntu` user), which files copied in as
`ubuntu` are.

### Rebuilds

A rebuild runs in the background, so you can log out; the site keeps serving the old
build until the new one is ready, then switches and restarts by itself. On the dev
instance (m6i.xlarge, builds on 3 vCPUs) an SGD rebuild takes about 20 minutes: ~7 for
QC, ~10 to precompute every gene. `./spellctl status sgd` shows progress, and
`/srv/spell/sgd/rebuild.log` keeps the full output. Errors stop the rebuild before anything
changes; only one rebuild per MOD runs at a time. Use `./spellctl rebuild sgd --wait` to
stay attached instead. Per-dataset results are in
`/srv/spell/sgd/validation_report.tsv`; QC results in `current/qc_report.tsv`.

## Admin tasks

| Task | Command |
| --- | --- |
| First set-up of a MOD | `./spellctl init wb ../configs/wb.example.toml`, fill `raw/`, `./spellctl rebuild wb` |
| Rebuild the image after a code update | `git pull`, `./spellctl image`, `./spellctl start sgd` |
| Start / stop an app | `./spellctl start sgd`, `./spellctl stop sgd` |

`start` recreates the container on the new image; `restart` would keep running the old one.

### Weekly gene name and GO update

Gene names, aliases and descriptions come from `SGD_features.tab`, and GO enrichment from
`go-basic.obo` and the SGD GAF. The `[updates]` section of `raw/mod.toml` lists their download
URLs (see `configs/sgd.toml`). `./spellctl update sgd` fetches them and, only if a file changed,
rebuilds (about 20 minutes; the site stays up). A download is swapped in only if it passes
checks (not empty, not an HTML error page, a complete gzip file, at least half the size of the
current file); the replaced file is kept as `<name>.prev`. If a rebuild is already running,
that week's update is skipped.

cron runs it weekly on the instance (`crontab -l`):

```
# spell-next: weekly gene name and GO update for SGD (Sundays 10:00 UTC = 2-3 am Pacific)
0 10 * * 0  /home/ubuntu/spell-next/deploy/spellctl update sgd >> /srv/spell/sgd/update.log 2>&1
```

`./spellctl status sgd` shows the last update; `/srv/spell/sgd/update.log` has the details.
To update by hand (e.g. after an urgent gene rename), run `./spellctl update sgd`. Files a
MOD edits by hand (expression files, `datasets.tsv`) are never downloaded.

### Code updates

`~/spell-next` on the instance is a git checkout of `main`. It pulls with a read-only
deploy key (`~/.ssh/spell_next_deploy`, used for `github.com` via `~/.ssh/config`; listed
under the repository's Settings → Deploy keys as "spell-dev (read-only)"), so nobody's
personal GitHub credentials are needed on the instance. To deploy pushed changes:

```bash
cd ~/spell-next && git pull && cd deploy
./spellctl image && ./spellctl start sgd
```

Data rebuilds are only needed when a change affects the build (search, QC or
precompute code); front-end and API changes need just the steps above. `deploy/.env`
holds instance-only settings (`SPELL_BIND=0.0.0.0` so the load balancer can reach the
app); it is not in git and `git pull` leaves it alone.

Containers restart automatically after a reboot (`restart: unless-stopped`). Each app
is limited to 2 CPUs and 6 GB of memory; builds get 3 CPUs and 10 GB
(`SPELL_BUILD_CPUS`, `SPELL_BUILD_MEM`).
Back up `/srv/spell/*/raw` and `/srv/spell/*/logs`; builds can be recreated from `raw/`.

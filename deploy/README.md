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

Ports are not public. With the instance's SSH key:

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

**Adding a dataset:** copy the PCL file to `/srv/spell/sgd/raw/expression/`, add
`<ID>\t<file name>` to `raw/dataset_list.txt` (validate prints the next free ID), add a
row to `raw/meta/datasets.tsv` (file name, GEO ID, PubMed ID, title, ... tags), then
`./spellctl validate sgd` and `./spellctl rebuild sgd`.

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
| Rebuild the image after a code update | `./spellctl image`, then `./spellctl restart sgd` |
| Start / stop an app | `./spellctl start sgd`, `./spellctl stop sgd` |

Containers restart automatically after a reboot (`restart: unless-stopped`). Each app
is limited to 2 CPUs and 6 GB of memory; builds get 3 CPUs and 10 GB
(`SPELL_BUILD_CPUS`, `SPELL_BUILD_MEM`).
Back up `/srv/spell/*/raw` and `/srv/spell/*/logs`; builds can be recreated from `raw/`.

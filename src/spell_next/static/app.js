'use strict';

const app = document.getElementById('app');
const tooltip = document.getElementById('tooltip');
const HEATMAP_DATASETS = 8;
const HEATMAP_GENES = 15;
const GO_GENES = 50;      // GO enrichment uses the top hits only
const CLIP = 2;   // heatmap colour saturates at +/- 2 standard deviations from the gene's mean
const CELL_H = 11;
const GAP = 4;    // space between the query genes and their hits

// Heatmap display options (like SPELL's "Additional Display Options"); remembered per browser.
const SCALES = {
  z: 'Standardized (z-score)',
  fold: 'Log change from gene mean',
  rank: 'Rank order',
};
const SCHEMES = {
  bluered: { name: 'Blue / red', low: [33, 102, 172], mid: [247, 247, 247], high: [178, 24, 43] },
  redgreen: { name: 'Red / green', low: [0, 170, 0], mid: [0, 0, 0], high: [220, 0, 0] },
  yellowblue: { name: 'Yellow / blue', low: [0, 60, 220], mid: [0, 0, 0], high: [235, 220, 0] },
};
const display = (() => {
  const saved = {};
  try { Object.assign(saved, JSON.parse(localStorage.getItem('spellDisplay') || '{}')); } catch (e) { /* ignore */ }
  return { scale: SCALES[saved.scale] ? saved.scale : 'z', scheme: SCHEMES[saved.scheme] ? saved.scheme : 'bluered' };
})();
function saveDisplay() {
  try { localStorage.setItem('spellDisplay', JSON.stringify(display)); } catch (e) { /* ignore */ }
}

const esc = (s) => String(s ?? '').replace(/[&<>"']/g,
  (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
const label = (g) => g.standard || g.systematic;
const sgdLink = (g) => `https://www.yeastgenome.org/locus/${encodeURIComponent(g.sgdid || g.systematic)}`;

async function api(path, options) {
  const res = await fetch(path, options);
  const body = await res.json().catch(() => ({}));
  if (!res.ok) {
    const detail = body.detail;
    throw new Error(typeof detail === 'string' ? detail : (detail && detail.message) || res.statusText);
  }
  return body;
}

// AG Grid: shared look and 10-row pages for every results table.
agGrid.ModuleRegistry.registerModules([agGrid.AllCommunityModule]);
const GRID_THEME = agGrid.themeQuartz.withPart(agGrid.colorSchemeVariable).withParams({ fontSize: 14, headerFontSize: 12 });
const darkQuery = window.matchMedia('(prefers-color-scheme: dark)');
const setGridMode = () => { document.body.dataset.agThemeMode = darkQuery.matches ? 'dark' : 'light'; };
setGridMode();
darkQuery.addEventListener('change', setGridMode);

function makeGrid(container, columnDefs, rowData) {
  return agGrid.createGrid(container, {
    theme: GRID_THEME,
    columnDefs,
    rowData,
    defaultColDef: { sortable: true, resizable: true, suppressMovable: true },
    pagination: true,
    paginationPageSize: 10,
    paginationPageSizeSelector: [10, 25, 50],
    domLayout: 'autoHeight',
    rowHeight: 60,
    suppressCellFocus: true,
    enableCellTextSelection: true,
  });
}

// "Download CSV" button for a grid; exports every row with raw (unformatted) values.
function csvButton(section, gridApi, fileName, columnKeys) {
  const btn = document.createElement('button');
  btn.type = 'button';
  btn.className = 'csv';
  btn.textContent = 'Download CSV';
  btn.addEventListener('click', () => gridApi.exportDataAsCsv({
    fileName,
    columnKeys,
    allColumns: true,
    processCellCallback: (p) => (Array.isArray(p.value)
      ? p.value.map((v) => (typeof v === 'object' ? `${v.term_id} ${v.name}` : v)).join('; ')
      : p.value),
  }));
  section.querySelector('h2').appendChild(btn);
}

const fileSlug = (genes) => genes.map(label).join('_').replace(/[^A-Za-z0-9_-]+/g, '').slice(0, 60) || 'query';

function cell(html) {
  const div = document.createElement('div');
  div.className = 'cell';
  div.innerHTML = `<div class="cell-inner">${html}</div>`;
  return div;
}

// "First author year · GEO · N conditions" without stray separators when parts are missing.
const datasetMeta = (d, extra = []) => [[d.first_author, d.year].filter(Boolean).join(' '), d.geo_id, ...extra]
  .filter(Boolean).map(esc).join(' · ');

function bar(text, fraction) {
  return cell(`<div class="bar"><span>${esc(text)}</span>
    <i style="width:${(100 * Math.max(0, Math.min(1, fraction))).toFixed(1)}%"></i></div>`);
}

function route() {
  const [path, qs] = location.hash.slice(1).split('?');
  const [, view, arg] = path.split('/').map((p) => decodeURIComponent(p || ''));
  const tags = (new URLSearchParams(qs || '').get('tags') || '').split('|').filter(Boolean);
  if (view === 'search' && arg) return showSearch(arg, tags);
  if (view === 'dataset' && arg) return showDataset(Number(arg));
  if (view === 'datasets') return showDatasets(tags);
  if (view === 'expression' && arg) return showListing(arg, tags);
  if (view === 'help') return showHelp();
  if (view === 'about') return showAbout();
  return showHome();
}

const searchHash = (query, tags = []) =>
  `#/search/${encodeURIComponent(query)}${tags.length ? `?tags=${encodeURIComponent(tags.join('|'))}` : ''}`;

let tagCache = null;
async function allTags() {
  if (!tagCache) tagCache = await api('/api/tags');
  return tagCache;
}

// Checkbox list of dataset tags; calls onChange(selected Set) when a box is toggled.
const TAGS_SHOWN = 15;

function tagPicker(container, tags, selected, onChange) {
  container.classList.add('tag-picker');
  const hidden = (i) => i >= TAGS_SHOWN && !selected.has(tags[i].tag);
  const extra = tags.filter((t, i) => hidden(i)).length;
  container.innerHTML = tags.map((t, i) => `<label class="tag${hidden(i) ? ' more' : ''}">
    <input type="checkbox" data-i="${i}" ${selected.has(t.tag) ? 'checked' : ''}>
    ${esc(t.tag)} <span class="sub">${t.count}</span></label>`).join('') +
    (extra ? `<button type="button" class="tag-toggle secondary" aria-expanded="false">Show all ${tags.length}
      tags</button>` : '');
  const toggle = container.querySelector('.tag-toggle');
  if (toggle) {
    toggle.addEventListener('click', () => {
      const open = container.classList.toggle('show-all');
      toggle.setAttribute('aria-expanded', open);
      toggle.textContent = open ? 'Show fewer tags' : `Show all ${tags.length} tags`;
    });
  }
  container.querySelectorAll('input').forEach((box) => box.addEventListener('change', () => {
    const tag = tags[Number(box.dataset.i)].tag;
    if (box.checked) selected.add(tag); else selected.delete(tag);
    onChange(selected);
  }));
}

function showHome() {
  document.getElementById('search-input').value = '';
  const examples = ['ACT1', 'GAL1', 'HSP104', 'RPL3 RPL4A RPS2 RPS3', 'CDC28 CLN1 CLN2 CLB5 CLB6'];
  app.innerHTML = `
    <h1>Find genes co-expressed with your query</h1>
    <p class="sub">Enter one gene for precomputed results, or several genes to search live for genes
      that share their expression pattern across about 800 yeast expression datasets.</p>
    <p>Examples: ${examples.map((e) => `<a href="#/search/${encodeURIComponent(e)}">${esc(e)}</a>`).join(' · ')}</p>
    <p class="sub">To see the expression of genes without a co-expression search, use a listing, e.g.
      <a href="#/expression/${encodeURIComponent('GAL1 GAL7 GAL10')}">GAL1 GAL7 GAL10</a>.
      Browse or filter all datasets by topic on the <a href="#/datasets">Datasets</a> page.</p>
    <p class="sub">Prototype rebuilt from SGD SPELL. Scores use corrected correlation arithmetic, so they
      can differ from spell.yeastgenome.org.</p>`;
}

async function showSearch(query, tags = []) {
  document.getElementById('search-input').value = query;
  app.innerHTML = `<p class="loading">Searching for ${esc(query)}…</p>`;
  let result;
  try {
    const names = query.split(/[\s,;]+/).filter(Boolean);
    const lookup = await api(`/api/genes/lookup?q=${encodeURIComponent(names.join(' '))}`);
    if (tags.length) {
      result = await api('/api/search', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ genes: names, limit: 100, datasets: 100, tags }),
      });
    } else if (lookup.genes.length === 1) {
      result = await api(`/api/genes/${lookup.genes[0].id}/coexpression?limit=100&datasets=100`);
      result.unknown = lookup.unknown;
      result.ambiguous = lookup.ambiguous;
    } else {
      result = await api('/api/search', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ genes: names, limit: 100, datasets: 100 }),
      });
    }
  } catch (err) {
    app.innerHTML = `<p class="notice">${esc(err.message)}</p>`;
    return;
  }
  renderResults(result, query);
}

function renderResults(r, query) {
  const notices = [];
  if (r.unknown && r.unknown.length) notices.push(`Not found: ${r.unknown.map(esc).join(', ')}`);
  for (const a of r.ambiguous || []) {
    notices.push(`${esc(a.name)} matches several genes: ` +
      a.candidates.map((c) => `<a href="#/search/${encodeURIComponent(c.systematic)}">${esc(label(c))}
        (${esc(c.systematic)})</a>`).join(', '));
  }
  const maxScore = Math.max(...r.genes.map((g) => g.score), 1e-9);
  const maxWeight = Math.max(...r.datasets.map((d) => d.weight), 1e-9);
  const single = r.query.length === 1;

  app.innerHTML = `
    <h1>${single ? 'Genes co-expressed with' : 'Genes co-expressed with this set'}</h1>
    <div class="chips">${r.query.map((g) => `<span class="chip" title="${esc(g.description)}">
      ${esc(label(g))}</span>`).join('')}</div>
    ${notices.map((n) => `<p class="notice">${n}</p>`).join('')}
    <p class="sub">${r.source === 'live' ? 'Live search' : 'Precomputed result'} ·
      top ${r.genes.length} genes · ${r.datasets.length} most informative datasets shown</p>
    ${r.tags && r.tags.length ? `<p class="notice">Searched only the ${r.datasets_searched} datasets tagged
      ${r.tags.map((t) => `<b>${esc(t)}</b>`).join(' or ')}.
      <a href="${searchHash(query)}">Search all datasets</a></p>` : ''}
    <p class="sub"><a href="#/expression/${encodeURIComponent(query)}${r.tags && r.tags.length ?
      `?tags=${encodeURIComponent(r.tags.join('|'))}` : ''}">Show the expression of the query genes only</a>
      (no co-expression search)</p>
    <details class="tagfilter">
      <summary>Restrict the search to datasets with these tags${r.tags && r.tags.length ?
        ` (${r.tags.length} selected)` : ''}</summary>
      <div id="tag-panel" class="tag-picker"><p class="loading">Loading tags…</p></div>
      <div class="tag-actions"><button type="button" id="tags-apply">Search selected tags</button>
        <button type="button" id="tags-clear" class="secondary">Clear</button></div>
    </details>
    <div class="grid">
      <section>
        <h2>Co-expressed genes</h2>
        <div id="genes-grid"></div>
      </section>
      <section>
        <h2>Most informative datasets</h2>
        <div id="datasets-grid"></div>
      </section>
    </div>
    <h2>GO enrichment of the top ${Math.min(GO_GENES, r.genes.length)} genes</h2>
    <div class="tabs" role="tablist">
      ${[['P', 'Process'], ['F', 'Function'], ['C', 'Component']].map(([a, l], i) =>
        `<button type="button" role="tab" data-aspect="${a}" aria-selected="${i === 0}">${l}</button>`).join('')}
    </div>
    <div id="go-results"><p class="loading">Calculating enrichment…</p></div>
    <h2>Expression in the most informative datasets</h2>
    <div id="heatmaps"></div>`;

  const slug = fileSlug(r.query);
  const genesGrid = makeGrid(document.getElementById('genes-grid'), [
    { headerName: '#', field: 'rank', width: 64, type: 'numericColumn' },
    { headerName: 'Systematic name', field: 'systematic', hide: true },
    { headerName: 'SGD ID', field: 'sgdid', hide: true },
    { headerName: 'Gene', field: 'label', width: 150, cellRenderer: (c) => cell(
      `<a href="#/search/${encodeURIComponent(c.data.systematic)}">${esc(c.value)}</a><br>
       <a class="sub nowrap" href="${sgdLink(c.data)}" target="_blank" rel="noopener">${esc(c.data.systematic)}&nbsp;↗</a>`) },
    { headerName: 'Score', field: 'score', width: 130, sort: 'desc',
      cellRenderer: (c) => bar(c.value.toFixed(2), c.value / maxScore) },
    { headerName: 'Description', field: 'description', flex: 1, minWidth: 220,
      cellRenderer: (c) => cell(`<div class="clamp" title="${esc(c.value)}">${esc(c.value)}</div>`) },
  ], r.genes.map((g, i) => ({ ...g, rank: i + 1, label: label(g) })));
  csvButton(document.getElementById('genes-grid').parentElement, genesGrid, `spell-next_${slug}_genes.csv`,
    ['rank', 'label', 'systematic', 'sgdid', 'score', 'description']);

  const datasetsGrid = makeGrid(document.getElementById('datasets-grid'), [
    { headerName: '#', field: 'rank', width: 64, type: 'numericColumn' },
    { headerName: 'Dataset ID', field: 'id', hide: true },
    { headerName: 'First author', field: 'first_author', hide: true },
    { headerName: 'Year', field: 'year', hide: true },
    { headerName: 'GEO', field: 'geo_id', hide: true },
    { headerName: 'PubMed', field: 'pubmed_id', hide: true },
    { headerName: 'Conditions', field: 'num_conds', hide: true },
    { headerName: 'Dataset', field: 'title', flex: 1, minWidth: 240, cellRenderer: (c) => cell(
      `<div class="clamp1"><a href="#/dataset/${c.data.id}" title="${esc(c.value)}">${esc(c.value)}</a></div>
       <span class="sub">${datasetMeta(c.data, [`${c.data.num_conds} conditions`])}</span>`) },
    { headerName: 'Weight', field: 'weight', width: 120, sort: 'desc',
      cellRenderer: (c) => bar(`${c.value.toFixed(1)}%`, c.value / maxWeight) },
  ], r.datasets.map((d, i) => ({ ...d, rank: i + 1 })));
  csvButton(document.getElementById('datasets-grid').parentElement, datasetsGrid, `spell-next_${slug}_datasets.csv`,
    ['rank', 'id', 'title', 'first_author', 'year', 'geo_id', 'pubmed_id', 'num_conds', 'weight']);

  const hitIds = r.genes.slice(0, GO_GENES).map((g) => g.id);
  document.querySelectorAll('.tabs button').forEach((b) => b.addEventListener('click', () => {
    document.querySelectorAll('.tabs button').forEach((x) => x.setAttribute('aria-selected', x === b));
    showEnrichment(hitIds, b.dataset.aspect, r.query);
  }));
  showEnrichment(hitIds, 'P', r.query);

  const picked = new Set(r.tags || []);
  allTags().then((tags) => tagPicker(document.getElementById('tag-panel'), tags, picked, () => {}))
    .catch((err) => { document.getElementById('tag-panel').innerHTML = `<p class="notice">${esc(err.message)}</p>`; });
  document.getElementById('tags-apply').addEventListener('click', () => {
    location.hash = searchHash(query, [...picked]);
  });
  document.getElementById('tags-clear').addEventListener('click', () => { location.hash = searchHash(query); });

  const genes = [...r.query, ...r.genes.slice(0, HEATMAP_GENES)];
  heatmapPager(document.getElementById('heatmaps'), r.datasets, genes, r.query.length);
}

// Heatmaps for a list of datasets, HEATMAP_DATASETS per page.
function heatmapPager(box, datasets, genes, nQuery) {
  const pages = Math.max(1, Math.ceil(datasets.length / HEATMAP_DATASETS));
  let page = 0;
  box.innerHTML = `
    <div class="display-options">
      <label>Scale <select data-opt="scale">${Object.entries(SCALES).map(([k, v]) =>
        `<option value="${k}" ${display.scale === k ? 'selected' : ''}>${esc(v)}</option>`).join('')}</select></label>
      <label>Colours <select data-opt="scheme">${Object.entries(SCHEMES).map(([k, v]) =>
        `<option value="${k}" ${display.scheme === k ? 'selected' : ''}>${esc(v.name)}</option>`).join('')}</select></label>
    </div>
    <div class="legend"></div>
    <div class="pager">
      <button type="button" data-step="-1">‹ Previous</button>
      <span class="pager-label"></span>
      <button type="button" data-step="1">Next ›</button>
    </div>
    <div class="heatmaps"></div>`;
  const grid = box.querySelector('.heatmaps');
  const legend = () => {
    const what = {
      z: `higher than the gene's average in that dataset (±${CLIP} standard deviations)`,
      fold: `higher than the gene's average in that dataset (±${CLIP} log units)`,
      rank: 'higher-ranked condition for that gene in that dataset',
    }[display.scale];
    box.querySelector('.legend').innerHTML = `<span>lower</span>
      <span class="ramp" style="background:${rampCss()}"></span>
      <span>${esc(what)}. Click a row to see that gene's full profile.</span>`;
  };
  box.querySelectorAll('.display-options select').forEach((sel) => sel.addEventListener('change', () => {
    display[sel.dataset.opt] = sel.value;
    saveDisplay();
    legend();
    show();
  }));
  legend();
  const show = () => {
    const first = page * HEATMAP_DATASETS;
    const shown = datasets.slice(first, first + HEATMAP_DATASETS);
    box.querySelector('.pager-label').textContent =
      `Datasets ${first + 1}–${first + shown.length} of ${datasets.length}`;
    box.querySelector('[data-step="-1"]').disabled = page === 0;
    box.querySelector('[data-step="1"]').disabled = page >= pages - 1;
    grid.innerHTML = '';
    shown.forEach((d) => drawHeatmap(grid, d, genes, nQuery));
  };
  box.querySelectorAll('.pager button').forEach((b) => b.addEventListener('click', () => {
    page = Math.min(pages - 1, Math.max(0, page + Number(b.dataset.step)));
    show();
  }));
  show();
}

async function showEnrichment(geneIds, aspect, queryGenes) {
  const box = document.getElementById('go-results');
  box.innerHTML = '<p class="loading">Calculating enrichment…</p>';
  let e;
  try {
    e = await api('/api/go-enrichment', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ gene_ids: geneIds, aspect, max_results: 100 }),
    });
  } catch (err) {
    box.innerHTML = `<p class="notice">${esc(err.message)}</p>`;
    return;
  }
  const sig = e.terms.filter((t) => t.fdr < 0.05);
  const p = (v) => (v < 1e-3 ? v.toExponential(1) : v.toFixed(3));
  const source = `GO ${esc(e.versions.go_ontology || '')}, SGD annotations ${esc(e.versions.go_annotations || '')};
    background: ${e.background_size} annotated genes`;
  if (!sig.length) {
    box.innerHTML = `<p class="sub">No ${esc(e.aspect.replace('_', ' '))} terms enriched at FDR &lt; 5%. ${source}</p>`;
    return;
  }
  const quickGO = (id) => `https://www.ebi.ac.uk/QuickGO/term/${encodeURIComponent(id)}`;
  box.innerHTML = `<h2 class="csv-row"></h2><div id="go-grid"></div>
    <p class="sub">Hypergeometric test, Bonferroni-corrected p and Benjamini-Hochberg FDR. Terms with the same
      genes from the list as a higher-ranked term are listed under it ("also"). ${source}</p>`;
  const goGrid = makeGrid(document.getElementById('go-grid'), [
    { headerName: 'GO ID', field: 'term_id', hide: true },
    { headerName: 'p', field: 'p_value', hide: true },
    { headerName: 'Also (same list genes)', field: 'also', hide: true },
    { headerName: 'GO term', field: 'name', flex: 2, minWidth: 260, cellRenderer: (c) => cell(
      `<a href="${quickGO(c.data.term_id)}" target="_blank" rel="noopener">${esc(c.value)}</a>
       <span class="sub nowrap">${esc(c.data.term_id)}</span>
       ${c.data.also.length ? `<div class="also clamp1" title="${esc(c.data.also.map((a) => a.name).join('; '))}">also:
         ${c.data.also.map((a) => `<a href="${quickGO(a.term_id)}" target="_blank" rel="noopener">${esc(a.name)}</a>`)
           .join(', ')}</div>` : ''}`) },
    { headerName: 'In list', field: 'list_count', width: 100, type: 'numericColumn',
      valueFormatter: (c) => `${c.value}/${c.data.list_size}` },
    { headerName: 'Genome', field: 'background_count', width: 110, type: 'numericColumn',
      valueFormatter: (c) => `${c.value}/${c.data.background_size}` },
    { headerName: 'Corrected p', field: 'p_bonferroni', width: 125, type: 'numericColumn',
      valueFormatter: (c) => p(c.value) },
    { headerName: 'FDR', field: 'fdr', width: 100, type: 'numericColumn', sort: 'asc',
      valueFormatter: (c) => p(c.value) },
    { headerName: 'Genes', field: 'genes', flex: 2, minWidth: 220, sortable: false,
      cellRenderer: (c) => cell(`<div class="clamp" title="${esc(c.value.join(', '))}">${c.value.map(esc).join(', ')}</div>`) },
  ], sig);
  csvButton(box, goGrid, `spell-next_${fileSlug(queryGenes)}_go_${aspect}.csv`,
    ['term_id', 'name', 'list_count', 'background_count', 'p_value', 'p_bonferroni', 'fdr', 'genes', 'also']);
}

async function drawHeatmap(box, dataset, genes, nQuery) {
  const card = document.createElement('div');
  card.className = 'heatmap';
  card.innerHTML = `<h3><a href="#/dataset/${dataset.id}">${esc(dataset.title)}</a></h3>
    <div class="meta">${esc(dataset.first_author || '')} ${esc(dataset.year || '')}</div>`;
  box.appendChild(card);
  let data;
  try {
    data = await api(`/api/datasets/${dataset.id}/expression?genes=${genes.map((g) => g.id).join(',')}`);
  } catch (err) {
    card.insertAdjacentHTML('beforeend', `<p class="notice">${esc(err.message)}</p>`);
    return;
  }
  const nCond = data.conditions.length;
  const z = data.rows.map((row) => zscores(row.values));
  const mapped = data.rows.map((row) => mapValues(row.values));
  const rowTop = (gi) => gi * CELL_H + (gi >= nQuery ? GAP : 0);

  // One canvas pixel per condition, stretched horizontally; rows keep CELL_H px so labels line up.
  const canvas = document.createElement('canvas');
  canvas.width = nCond;
  canvas.height = rowTop(genes.length);
  canvas.style.height = `${canvas.height}px`;
  canvas.setAttribute('role', 'img');
  canvas.setAttribute('aria-label', `Expression heatmap for ${dataset.title}`);
  const ctx = canvas.getContext('2d');
  mapped.forEach((vals, gi) => {
    for (let c = 0; c < nCond; c++) {
      ctx.fillStyle = colour(vals ? vals[c] : null);
      ctx.fillRect(c, rowTop(gi), 1, CELL_H - 1);
    }
  });
  canvas.addEventListener('mousemove', (ev) => {
    const rect = canvas.getBoundingClientRect();
    const c = Math.floor((ev.clientX - rect.left) * nCond / rect.width);
    const y = ev.clientY - rect.top;
    const gi = y < nQuery * CELL_H ? Math.floor(y / CELL_H) : Math.floor((y - GAP) / CELL_H);
    const row = data.rows[gi];
    if (!row || c < 0 || c >= nCond) { tooltip.hidden = true; return; }
    const v = row.values ? row.values[c] : null;
    tooltip.textContent = `${label(genes[gi])} · ${data.conditions[c]} · ` +
      (v === null ? 'no data' : `${v.toFixed(2)} from mean (z ${z[gi][c].toFixed(1)})`) + ' · click for profile';
    tooltip.style.left = `${ev.clientX + 12}px`;
    tooltip.style.top = `${ev.clientY + 12}px`;
    tooltip.hidden = false;
  });
  canvas.addEventListener('mouseleave', () => { tooltip.hidden = true; });
  canvas.addEventListener('click', (ev) => {
    const y = ev.clientY - canvas.getBoundingClientRect().top;
    const gi = y < nQuery * CELL_H ? Math.floor(y / CELL_H) : Math.floor((y - GAP) / CELL_H);
    if (genes[gi]) openProfile(dataset, genes[gi]);
  });

  const labels = document.createElement('div');
  labels.className = 'labels';
  labels.innerHTML = genes.map((g, gi) => `<div style="top:${rowTop(gi)}px"
    class="${gi < nQuery ? 'q' : ''}">${esc(label(g))}</div>`).join('');
  labels.style.height = `${canvas.height}px`;
  const body = document.createElement('div');
  body.className = 'hm-body';
  body.append(labels, canvas);
  card.appendChild(body);
}

// Full profile of one gene in one dataset: a labelled bar per condition (replaces SPELL's popup).
async function openProfile(dataset, gene) {
  const dialog = document.getElementById('profile');
  const body = dialog.querySelector('.profile-body');
  dialog.querySelector('.profile-title').innerHTML = `${esc(label(gene))}
    <span class="sub">${esc(gene.systematic)}</span>`;
  dialog.querySelector('.profile-dataset').innerHTML =
    `<a href="#/dataset/${dataset.id}">${esc(dataset.title)}</a>`;
  body.innerHTML = '<p class="loading">Loading…</p>';
  if (!dialog.open) dialog.showModal();
  let data;
  try {
    data = await api(`/api/datasets/${dataset.id}/expression?genes=${gene.id}&raw=true`);
  } catch (err) {
    body.innerHTML = `<p class="notice">${esc(err.message)}</p>`;
    return;
  }
  const row = data.rows[0];
  if (!row.values) {
    body.innerHTML = '<p class="notice">This gene was not measured in this dataset.</p>';
    return;
  }
  body.innerHTML = profileChart(data.conditions, row.values, row.raw);
}

function profileChart(conditions, values, raw) {
  const W = 680, ROW = 16, LABEL_W = 260, VALUE_W = 64, TOP = 22;
  const plotW = W - LABEL_W - VALUE_W - 16;
  const zero = LABEL_W + 8 + plotW / 2;
  const maxAbs = Math.max(1, ...values.filter((v) => v !== null).map(Math.abs));
  const x = (v) => zero + (v / maxAbs) * (plotW / 2);
  const H = TOP + conditions.length * ROW + 8;
  const short = (t) => (t.length > 42 ? `${t.slice(0, 41)}…` : t);
  const rows = conditions.map((cond, i) => {
    const y = TOP + i * ROW;
    const v = values[i];
    const tip = `${cond}: ${v === null ? 'no data' : `${v.toFixed(2)} from mean, value ${raw[i].toFixed(2)}`}`;
    const barEl = v === null ? '' : `<rect x="${Math.min(zero, x(v))}" y="${y + 3}" width="${Math.abs(x(v) - zero)}"
      height="${ROW - 6}" fill="${colour(v >= 0 ? 1 : -1)}"/>`;
    return `<g><title>${esc(tip)}</title>
      <text x="${LABEL_W}" y="${y + 12}" text-anchor="end" class="lbl">${esc(short(cond))}</text>${barEl}
      <text x="${W - 4}" y="${y + 12}" text-anchor="end" class="val">${v === null ? '—' : v.toFixed(2)}</text></g>`;
  }).join('');
  return `<p class="sub">Difference from the gene's mean in this dataset (log scale); hover a row for the
      stored value. ${conditions.length} conditions.</p>
    <svg class="profile-chart" viewBox="0 0 ${W} ${H}" width="100%" role="img"
      aria-label="Expression profile across ${conditions.length} conditions">
      <text x="${zero}" y="12" text-anchor="middle" class="axis">0</text>
      <text x="${zero - plotW / 2}" y="12" class="axis">−${maxAbs.toFixed(1)}</text>
      <text x="${zero + plotW / 2}" y="12" text-anchor="end" class="axis">+${maxAbs.toFixed(1)}</text>
      <line x1="${zero}" x2="${zero}" y1="${TOP - 4}" y2="${H - 4}" class="zero"/>
      ${rows}
    </svg>`;
}

// Standardize one gene's (already mean-centred) values so every dataset uses the same colour scale.
function zscores(values) {
  if (!values) return null;
  const present = values.filter((v) => v !== null);
  if (present.length < 2) return values.map(() => null);
  const sd = Math.sqrt(present.reduce((s, v) => s + v * v, 0) / present.length);
  return values.map((v) => (v === null ? null : sd > 0 ? v / sd : 0));
}

// Map one gene's mean-centred values to [-1, 1] for colouring, per display.scale.
function mapValues(values) {
  if (!values) return null;
  const clip = (v) => Math.max(-1, Math.min(1, v));
  if (display.scale === 'fold') return values.map((v) => (v === null ? null : clip(v / CLIP)));
  if (display.scale === 'rank') {
    const order = values.map((v, i) => [v, i]).filter(([v]) => v !== null).sort((a, b) => a[0] - b[0]);
    const out = values.map(() => null);
    const n = order.length;
    order.forEach(([, i], r) => { out[i] = n > 1 ? (2 * r) / (n - 1) - 1 : 0; });
    return out;
  }
  const z = zscores(values);
  return z.map((v) => (v === null ? null : clip(v / CLIP)));
}

function colour(t) {
  if (t === null || t === undefined) return 'rgb(150,150,150)';
  const { low, mid, high } = SCHEMES[display.scheme];
  const end = t < 0 ? low : high;
  const f = Math.min(1, Math.abs(t));
  const mix = (k) => Math.round(mid[k] + (end[k] - mid[k]) * f);
  return `rgb(${mix(0)},${mix(1)},${mix(2)})`;
}

const rampCss = () => `linear-gradient(90deg, ${colour(-1)}, ${colour(0)}, ${colour(1)})`;

// Listing mode: the query genes' expression in every dataset that measures them (SPELL's "Show Expression Levels").
async function showListing(query, tags) {
  document.getElementById('search-input').value = query;
  app.innerHTML = `<p class="loading">Loading expression of ${esc(query)}…</p>`;
  let r;
  try {
    const qs = `genes=${encodeURIComponent(query)}${tags.length ? `&tags=${encodeURIComponent(tags.join('|'))}` : ''}`;
    r = await api(`/api/expression-listing?${qs}`);
  } catch (err) {
    app.innerHTML = `<p class="notice">${esc(err.message)}</p>`;
    return;
  }
  const n = r.query.length;
  app.innerHTML = `
    <h1>Expression of</h1>
    <div class="chips">${r.query.map((g) => `<span class="chip" title="${esc(g.description)}">
      ${esc(label(g))}</span>`).join('')}</div>
    ${r.unknown.length ? `<p class="notice">Not found: ${r.unknown.map(esc).join(', ')}</p>` : ''}
    <p class="sub">${r.datasets.length} datasets measure ${n > 1 ? 'these genes' : 'this gene'}${tags.length ?
      ` among those tagged ${tags.map(esc).join(' or ')}` : ''}. Datasets measuring the most query genes come
      first, then those where the genes vary most. No co-expression search is run;
      <a href="${searchHash(query, tags)}">search for co-expressed genes</a>.</p>
    <section><h2>Datasets</h2><div id="listing-grid"></div></section>
    <h2>Expression</h2>
    <div id="heatmaps"></div>`;
  const listingGrid = makeGrid(document.getElementById('listing-grid'), [
    { headerName: 'Dataset ID', field: 'id', hide: true },
    { headerName: 'GEO', field: 'geo_id', hide: true },
    { headerName: 'Dataset', field: 'title', flex: 1, minWidth: 260, cellRenderer: (c) => cell(
      `<div class="clamp1"><a href="#/dataset/${c.data.id}" title="${esc(c.value)}">${esc(c.value)}</a></div>
       <span class="sub">${datasetMeta(c.data)}</span>`) },
    { headerName: 'Genes measured', field: 'genes_measured', width: 150, type: 'numericColumn',
      valueFormatter: (c) => `${c.value}/${n}` },
    { headerName: 'Conditions', field: 'num_conds', width: 120, type: 'numericColumn' },
    { headerName: 'Mean variance', field: 'variance', width: 140, type: 'numericColumn',
      valueFormatter: (c) => c.value.toFixed(3) },
  ], r.datasets);
  csvButton(document.getElementById('listing-grid').parentElement, listingGrid,
    `spell-next_${fileSlug(r.query)}_listing.csv`, ['id', 'title', 'geo_id', 'genes_measured', 'num_conds', 'variance']);
  heatmapPager(document.getElementById('heatmaps'), r.datasets, r.query, n);
}

async function showHelp() {
  app.innerHTML = '<p class="loading">Loading…</p>';
  const i = await api('/api/info').catch(() => ({}));
  app.innerHTML = `<article class="doc">
    <h1>Help</h1>
    <h2>Searching</h2>
    <p>Type one or more gene names in the search box: standard names (ACT1), systematic names (YFL039C) or
      aliases, separated by spaces or commas (up to ${esc(i.max_query_genes || 50)} genes). If an alias matches
      several genes, you are asked to pick one.</p>
    <ul>
      <li><b>One gene</b>: results are precomputed and appear instantly.</li>
      <li><b>Several genes</b>: the search runs live and finds genes that share the expression pattern of the
        whole set. Pick genes you expect to act together; unrelated genes give weak results.</li>
    </ul>
    <h2>Reading the results</h2>
    <ul>
      <li><b>Score</b>: how strongly a gene is co-expressed with the query, summed over datasets (higher is
        stronger). Scores compare genes within one search, not across searches.</li>
      <li><b>Dataset weight</b>: how informative each dataset is for this query. For several genes it is how
        consistently the query genes are co-expressed in that dataset; for one gene, how much that gene varies.
        Weights add up to 100%.</li>
      <li><b>Heatmaps</b>: each row is a gene, each column a condition. Colours show expression relative to the
        gene's own average in that dataset (blue lower, red higher). Use Previous/Next to page through all
        listed datasets; click a row to see that gene's full profile with condition names.</li>
      <li><b>GO enrichment</b>: Gene Ontology terms over-represented among the top 50 genes (hypergeometric test;
        Bonferroni-corrected p and false discovery rate). Terms whose genes from the list are identical to a
        higher-ranked term are listed under it as "also".</li>
    </ul>
    <h2>Restricting by dataset tags</h2>
    <p>Datasets are tagged by topic (for example stress, heat shock, carbon utilization). On a results page open
      <i>Restrict the search to datasets with these tags</i>, or pick tags on the
      <a href="#/datasets">Datasets</a> page and search there. A dataset is used if it carries any selected tag.
      Unlike the original SPELL, which searched every dataset and only hid the others from the list, the search
      itself is limited to the tagged datasets, so the gene ranking reflects only them.</p>
    <h2>Expression without a search</h2>
    <p>A listing shows your genes' expression in every dataset that measures them, without ranking other genes:
      use <i>Show the expression of the query genes only</i> on a results page, or go to
      <code>#/expression/GAL1 GAL7 GAL10</code>. Tags restrict listings too.</p>
    <h2>Display options</h2>
    <p>Above the heatmaps, <i>Scale</i> switches between standardized values (z-scores), log change from the
      gene's mean, and rank order; <i>Colours</i> offers blue/red, red/green and yellow/blue. Your choice is
      remembered in this browser.</p>
    <h2>Downloading</h2>
    <p>Each table has a <i>Download CSV</i> button that saves every row with unrounded values. Dataset pages offer
      the whole expression matrix as searched (after any log transform) as a tab-separated file.</p>
    <h2>Embedding results on other sites</h2>
    <p>Other pages can show a compact result in an iframe:
      <code>&lt;iframe src="${esc(location.origin)}/embed?gene=ACT1&amp;genes=10&amp;sets=5"&gt;&lt;/iframe&gt;</code>
      (<code>genes</code>: co-expressed genes shown, <code>sets</code>: datasets shown).</p>
    <h2>Usage logging</h2>
    <p>Searches are logged by month (date, time, search type, query genes and tags) to report usage. No IP
      addresses or other personal data are recorded.</p>
    <h2>How results differ from the original SPELL</h2>
    <ul>
      <li>Correlations are computed in double precision. The original computed them in single precision with a
        formula that loses accuracy for some genes, which could change rankings.</li>
      <li>Data files are read correctly: empty cells no longer shift values into the wrong condition, and two
        datasets that were listed twice are counted once.</li>
      <li>Datasets stored without a log transform are log-transformed before searching
        (${esc((i.qc_scales || {}).ratio ?? '?')} of ratios, ${esc((i.qc_scales || {}).intensity ?? '?')} of
        intensities or counts). ${esc(i.qc_flagged ?? '?')} datasets are flagged for curator review on their
        dataset pages.</li>
      <li>GO enrichment uses current GO and SGD annotations.</li>
    </ul>
  </article>`;
}

async function showAbout() {
  app.innerHTML = '<p class="loading">Loading…</p>';
  const i = await api('/api/info').catch(() => ({ versions: {} }));
  const v = i.versions || {};
  const titleSearch = 'https://pubmed.ncbi.nlm.nih.gov/?term=' + encodeURIComponent(
    'Exploring the functional landscape of gene expression: directed search of large microarray compendia');
  app.innerHTML = `<article class="doc">
    <h1>About SPELL next</h1>
    <p>SPELL (Serial Pattern of Expression Levels Locator) is a query-driven search engine for large gene
      expression compendia. Given a set of query genes, it weights each dataset by how well the query genes are
      co-expressed in it, then ranks every other gene by its weighted correlation with the query. Negative
      correlations are treated as no correlation. This prototype re-implements SGD's SPELL in Python with a new
      web interface, as a candidate for hosting within the Alliance of Genome Resources.</p>
    <h2>Citation</h2>
    <p>The SPELL algorithm is described in: Hibbs MA, et al. Exploring the functional landscape of gene
      expression: directed search of large microarray compendia. <i>Bioinformatics</i>, 2007.
      <a href="${titleSearch}" target="_blank" rel="noopener">Find on PubMed ↗</a></p>
    <h2>Data</h2>
    <dl class="meta">
      <dt>Expression datasets</dt><dd>${esc(i.datasets ?? '—')} (${esc(i.conditions ?? '—')} conditions), from
        SGD SPELL; ${esc(i.datasets_with_metadata ?? '—')} with publication details</dd>
      <dt>Genes</dt><dd>${esc(i.genes ?? '—')}; names, aliases and descriptions from SGD</dd>
      <dt>Gene Ontology</dt><dd>${esc(v.go_ontology || '—')}</dd>
      <dt>GO annotations</dt><dd>SGD, ${esc(v.go_annotations || '—')}</dd>
    </dl>
    <h2>Credits</h2>
    <p>The original SPELL was developed by the Troyanskaya lab at Princeton University and the Hibbs lab at The
      Jackson Laboratory, and run at SGD as spell.yeastgenome.org. Questions about this prototype:
      <a href="mailto:sgd-helpdesk@lists.stanford.edu">sgd-helpdesk@lists.stanford.edu</a>.</p>
  </article>`;
}

async function showDatasets(initialTags) {
  app.innerHTML = '<p class="loading">Loading datasets…</p>';
  let datasets, tags;
  try {
    [datasets, tags] = await Promise.all([api('/api/datasets'), allTags()]);
  } catch (err) {
    app.innerHTML = `<p class="notice">${esc(err.message)}</p>`;
    return;
  }
  const selected = new Set(initialTags);
  app.innerHTML = `
    <h1>Datasets</h1>
    <p class="sub">${datasets.length} expression datasets. Pick tags to list matching datasets (a dataset
      matches if it carries any selected tag), then search genes within them.</p>
    <div class="tag-picker" id="ds-tags"></div>
    <div class="ds-tools">
      <input type="search" id="ds-filter" placeholder="Filter by title, author, GEO ID…" aria-label="Filter datasets">
      <form id="ds-search" class="search">
        <input type="search" id="ds-genes" placeholder="Gene name(s) to search in the selected datasets"
          aria-label="Genes to search in the selected datasets">
        <button type="submit">Search</button>
      </form>
    </div>
    <p class="sub" id="ds-count"></p>
    <div id="ds-grid"></div>`;
  const grid = agGrid.createGrid(document.getElementById('ds-grid'), {
    theme: GRID_THEME,
    rowData: datasets,
    columnDefs: [
      { headerName: 'Dataset', field: 'title', flex: 2, minWidth: 280, cellRenderer: (c) => cell(
        `<div class="clamp1"><a href="#/dataset/${c.data.id}" title="${esc(c.value)}">${esc(c.value)}</a></div>
         <span class="sub">${esc([c.data.first_author, c.data.year, c.data.journal].filter(Boolean).join(' · '))}</span>`) },
      { headerName: 'GEO', field: 'geo_id', width: 140 },
      { headerName: 'Year', field: 'year', width: 90, type: 'numericColumn' },
      { headerName: 'Conditions', field: 'num_conds', width: 120, type: 'numericColumn' },
      { headerName: 'Genes', field: 'num_genes', width: 100, type: 'numericColumn' },
      { headerName: 'Tags', field: 'tags', flex: 1, minWidth: 180, sortable: false,
        valueFormatter: (c) => c.value.join(', '),
        cellRenderer: (c) => cell(`<div class="clamp">${c.value.map(esc).join(', ') || '—'}</div>`) },
      { headerName: 'QC', field: 'warnings', width: 80, sortable: true,
        valueFormatter: (c) => (c.value ? 'review' : ''),
        cellRenderer: (c) => cell(c.value ? `<span class="flag" title="${esc(c.value)}">review</span>` : '') },
    ],
    defaultColDef: { sortable: true, resizable: true, suppressMovable: true },
    pagination: true,
    paginationPageSize: 10,
    paginationPageSizeSelector: [10, 25, 50, 100],
    domLayout: 'autoHeight',
    rowHeight: 60,
    suppressCellFocus: true,
    enableCellTextSelection: true,
    isExternalFilterPresent: () => selected.size > 0,
    doesExternalFilterPass: (node) => node.data.tags.some((t) => selected.has(t)),
    onModelUpdated: (e) => {
      document.getElementById('ds-count').textContent =
        `${e.api.getDisplayedRowCount()} of ${datasets.length} datasets shown`;
    },
  });
  tagPicker(document.getElementById('ds-tags'), tags, selected, () => {
    grid.onFilterChanged();
    history.replaceState(null, '', selected.size ? `#/datasets?tags=${encodeURIComponent([...selected].join('|'))}`
      : '#/datasets');
  });
  document.getElementById('ds-filter').addEventListener('input', (e) => {
    grid.setGridOption('quickFilterText', e.target.value);
  });
  document.getElementById('ds-search').addEventListener('submit', (e) => {
    e.preventDefault();
    const q = document.getElementById('ds-genes').value.trim();
    if (q) location.hash = searchHash(q, [...selected]);
  });
}

async function showDataset(id) {
  app.innerHTML = '<p class="loading">Loading dataset…</p>';
  let d;
  try {
    d = await api(`/api/datasets/${id}`);
  } catch (err) {
    app.innerHTML = `<p class="notice">${esc(err.message)}</p>`;
    return;
  }
  const links = [
    d.geo_id && `<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${esc(d.geo_id)}" target="_blank"
      rel="noopener">${esc(d.geo_id)} ↗</a>`,
    d.pubmed_id && `<a href="https://pubmed.ncbi.nlm.nih.gov/${d.pubmed_id}/" target="_blank"
      rel="noopener">PubMed ${d.pubmed_id} ↗</a>`,
  ].filter(Boolean).join(' · ');
  app.innerHTML = `
    <p><a href="javascript:history.back()">← Back</a></p>
    <h1>${esc(d.title)}</h1>
    <p class="sub">${esc(d.authors || d.first_author || '')}${d.journal ? ' · ' + esc(d.journal) : ''}
      ${d.year ? ' · ' + esc(d.year) : ''}</p>
    ${d.has_metadata ? '' : '<p class="notice">Publication details for this dataset are not yet loaded.</p>'}
    <dl class="meta">
      <dt>Links</dt><dd>${links || '—'}</dd>
      <dt>Conditions</dt><dd>${d.num_conds}</dd>
      <dt>Genes measured</dt><dd>${d.num_genes}</dd>
      <dt>Tags</dt><dd>${(d.tags || '').split('|').filter(Boolean).map((t) =>
        `<a href="#/datasets?tags=${encodeURIComponent(t)}">${esc(t)}</a>`).join(', ') || '—'}</dd>
      <dt>File</dt><dd>${esc(d.filename)}</dd>
      <dt>Download</dt><dd><a href="/api/datasets/${d.id}/download" download>Expression data as searched
        (tab-separated)</a></dd>
    </dl>
    ${d.qc ? qcSection(d.qc) : ''}
    ${d.description ? `<h2>Summary</h2><p>${esc(d.description)}</p>` : ''}
    <h2>Conditions</h2>
    <ol class="conds">${d.conditions.map((c) => `<li>${esc(c)}</li>`).join('')}</ol>`;
}

function qcSection(qc) {
  const split = (s) => (s || '').split('; ').filter(Boolean);
  const actions = split(qc.actions);
  const warnings = split(qc.warnings);
  const fmt = (v) => (v === null || v === undefined ? '—' : v.toFixed(2));
  return `<h2>Data quality</h2>
    ${warnings.map((w) => `<p class="notice">${esc(w)}</p>`).join('')}
    <dl class="meta">
      <dt>Value scale</dt><dd>${esc(qc.scale)}</dd>
      <dt>Processing</dt><dd>${actions.length ? actions.map(esc).join('; ') : 'none (values used as published)'}</dd>
      <dt>Condition shift</dt><dd>${fmt(qc.cond_median_spread)} (largest difference between condition medians)</dd>
      <dt>Baseline co-expression</dt><dd>${fmt(qc.z_mean_after)} mean Fisher z over all gene pairs
        ${qc.z_mean_before !== qc.z_mean_after ? `(${fmt(qc.z_mean_before)} before processing)` : ''}</dd>
    </dl>`;
}

document.getElementById('search-form').addEventListener('submit', (ev) => {
  ev.preventDefault();
  const q = document.getElementById('search-input').value.trim();
  if (q) location.hash = `#/search/${encodeURIComponent(q)}`;
});
window.addEventListener('hashchange', route);
route();

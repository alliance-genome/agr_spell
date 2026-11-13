'use strict';

const app = document.getElementById('app');
const tooltip = document.getElementById('tooltip');
const HEATMAP_DATASETS = 8;
const HEATMAP_GENES = 15;
const CLIP = 2;   // heatmap colour saturates at +/- 2 standard deviations from the gene's mean
const CELL_H = 11;
const GAP = 4;    // space between the query genes and their hits

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

function cell(html) {
  const div = document.createElement('div');
  div.className = 'cell';
  div.innerHTML = html;
  return div;
}

function bar(text, fraction) {
  return cell(`<div class="bar"><span>${esc(text)}</span>
    <i style="width:${(100 * Math.max(0, Math.min(1, fraction))).toFixed(1)}%"></i></div>`);
}

function route() {
  const hash = decodeURIComponent(location.hash.slice(1));
  const [, view, arg] = hash.split('/');
  if (view === 'search' && arg) return showSearch(arg);
  if (view === 'dataset' && arg) return showDataset(Number(arg));
  return showHome();
}

function showHome() {
  document.getElementById('search-input').value = '';
  const examples = ['ACT1', 'GAL1', 'HSP104', 'RPL3 RPL4A RPS2 RPS3', 'CDC28 CLN1 CLN2 CLB5 CLB6'];
  app.innerHTML = `
    <h1>Find genes co-expressed with your query</h1>
    <p class="sub">Enter one gene for precomputed results, or several genes to search live for genes
      that share their expression pattern across about 800 yeast expression datasets.</p>
    <p>Examples: ${examples.map((e) => `<a href="#/search/${encodeURIComponent(e)}">${esc(e)}</a>`).join(' · ')}</p>
    <p class="sub">Prototype rebuilt from SGD SPELL. Scores use corrected correlation arithmetic, so they
      can differ from spell.yeastgenome.org.</p>`;
}

async function showSearch(query) {
  document.getElementById('search-input').value = query;
  app.innerHTML = `<p class="loading">Searching for ${esc(query)}…</p>`;
  let result;
  try {
    const names = query.split(/[\s,;]+/).filter(Boolean);
    const lookup = await api(`/api/genes/lookup?q=${encodeURIComponent(names.join(' '))}`);
    if (lookup.genes.length === 1) {
      result = await api(`/api/genes/${lookup.genes[0].id}/coexpression?limit=50&datasets=20`);
      result.unknown = lookup.unknown;
      result.ambiguous = lookup.ambiguous;
    } else {
      result = await api('/api/search', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ genes: names, limit: 50, datasets: 20 }),
      });
    }
  } catch (err) {
    app.innerHTML = `<p class="notice">${esc(err.message)}</p>`;
    return;
  }
  renderResults(result);
}

function renderResults(r) {
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
    <h2>GO enrichment of the top ${r.genes.length} genes</h2>
    <div class="tabs" role="tablist">
      ${[['P', 'Process'], ['F', 'Function'], ['C', 'Component']].map(([a, l], i) =>
        `<button type="button" role="tab" data-aspect="${a}" aria-selected="${i === 0}">${l}</button>`).join('')}
    </div>
    <div id="go-results"><p class="loading">Calculating enrichment…</p></div>
    <h2>Expression in the top datasets</h2>
    <div class="legend"><span>lower</span><span class="ramp"></span>
      <span>higher than the gene's average in that dataset (±${CLIP} standard deviations)</span></div>
    <div id="heatmaps" class="heatmaps"></div>`;

  makeGrid(document.getElementById('genes-grid'), [
    { headerName: '#', field: 'rank', width: 64, type: 'numericColumn' },
    { headerName: 'Gene', field: 'label', width: 150, cellRenderer: (c) => cell(
      `<a href="#/search/${encodeURIComponent(c.data.systematic)}">${esc(c.value)}</a><br>
       <a class="sub nowrap" href="${sgdLink(c.data)}" target="_blank" rel="noopener">${esc(c.data.systematic)}&nbsp;↗</a>`) },
    { headerName: 'Score', field: 'score', width: 130, sort: 'desc',
      cellRenderer: (c) => bar(c.value.toFixed(2), c.value / maxScore) },
    { headerName: 'Description', field: 'description', flex: 1, minWidth: 220,
      cellRenderer: (c) => cell(`<div class="clamp" title="${esc(c.value)}">${esc(c.value)}</div>`) },
  ], r.genes.map((g, i) => ({ ...g, rank: i + 1, label: label(g) })));

  makeGrid(document.getElementById('datasets-grid'), [
    { headerName: '#', field: 'rank', width: 64, type: 'numericColumn' },
    { headerName: 'Dataset', field: 'title', flex: 1, minWidth: 240, cellRenderer: (c) => cell(
      `<div class="clamp1"><a href="#/dataset/${c.data.id}" title="${esc(c.value)}">${esc(c.value)}</a></div>
       <span class="sub">${esc([c.data.first_author, c.data.year].filter(Boolean).join(' '))}
       ${c.data.geo_id ? ' · ' + esc(c.data.geo_id) : ''} · ${c.data.num_conds} conditions</span>`) },
    { headerName: 'Weight', field: 'weight', width: 120, sort: 'desc',
      cellRenderer: (c) => bar(`${c.value.toFixed(1)}%`, c.value / maxWeight) },
  ], r.datasets.map((d, i) => ({ ...d, rank: i + 1 })));

  const hitIds = r.genes.map((g) => g.id);
  document.querySelectorAll('.tabs button').forEach((b) => b.addEventListener('click', () => {
    document.querySelectorAll('.tabs button').forEach((x) => x.setAttribute('aria-selected', x === b));
    showEnrichment(hitIds, b.dataset.aspect);
  }));
  showEnrichment(hitIds, 'P');

  const genes = [...r.query, ...r.genes.slice(0, HEATMAP_GENES)];
  const box = document.getElementById('heatmaps');
  r.datasets.slice(0, HEATMAP_DATASETS).forEach((d) => drawHeatmap(box, d, genes, r.query.length));
}

async function showEnrichment(geneIds, aspect) {
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
  box.innerHTML = `<div id="go-grid"></div>
    <p class="sub">Hypergeometric test, Bonferroni-corrected p and Benjamini-Hochberg FDR. Terms with the same
      genes from the list as a higher-ranked term are listed under it ("also"). ${source}</p>`;
  makeGrid(document.getElementById('go-grid'), [
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
  const rowTop = (gi) => gi * CELL_H + (gi >= nQuery ? GAP : 0);

  // One canvas pixel per condition, stretched horizontally; rows keep CELL_H px so labels line up.
  const canvas = document.createElement('canvas');
  canvas.width = nCond;
  canvas.height = rowTop(genes.length);
  canvas.style.height = `${canvas.height}px`;
  canvas.setAttribute('role', 'img');
  canvas.setAttribute('aria-label', `Expression heatmap for ${dataset.title}`);
  const ctx = canvas.getContext('2d');
  z.forEach((vals, gi) => {
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
      (v === null ? 'no data' : `${v.toFixed(2)} from mean (z ${z[gi][c].toFixed(1)})`);
    tooltip.style.left = `${ev.clientX + 12}px`;
    tooltip.style.top = `${ev.clientY + 12}px`;
    tooltip.hidden = false;
  });
  canvas.addEventListener('mouseleave', () => { tooltip.hidden = true; });

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

// Standardize one gene's (already mean-centred) values so every dataset uses the same colour scale.
function zscores(values) {
  if (!values) return null;
  const present = values.filter((v) => v !== null);
  if (present.length < 2) return values.map(() => null);
  const sd = Math.sqrt(present.reduce((s, v) => s + v * v, 0) / present.length);
  return values.map((v) => (v === null ? null : sd > 0 ? v / sd : 0));
}

function colour(v) {
  if (v === null || v === undefined) return 'rgb(150,150,150)';
  const t = Math.max(-1, Math.min(1, v / CLIP));
  const mix = (a, b, f) => Math.round(a + (b - a) * f);
  const white = [247, 247, 247];
  const end = t < 0 ? [33, 102, 172] : [178, 24, 43];
  const f = Math.abs(t);
  return `rgb(${mix(white[0], end[0], f)},${mix(white[1], end[1], f)},${mix(white[2], end[2], f)})`;
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
      <dt>Tags</dt><dd>${esc((d.tags || '').split('|').filter(Boolean).join(', ') || '—')}</dd>
      <dt>File</dt><dd>${esc(d.filename)}</dd>
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

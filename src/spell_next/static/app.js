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
        <div class="table-wrap"><table>
          <thead><tr><th>#</th><th>Gene</th><th>Score</th><th>Description</th></tr></thead>
          <tbody>${r.genes.map((g, i) => `
            <tr><td class="num">${i + 1}</td>
              <td><a href="#/search/${encodeURIComponent(g.systematic)}">${esc(label(g))}</a>
                <br><a class="sub nowrap" href="${sgdLink(g)}" target="_blank" rel="noopener">${esc(g.systematic)}&nbsp;↗</a></td>
              <td class="bar"><span>${g.score.toFixed(2)}</span>
                <i style="width:${(100 * g.score / maxScore).toFixed(1)}%"></i></td>
              <td class="desc"><div title="${esc(g.description)}">${esc(g.description)}</div></td></tr>`).join('')}
          </tbody></table></div>
      </section>
      <section>
        <h2>Most informative datasets</h2>
        <div class="table-wrap"><table>
          <thead><tr><th>#</th><th>Dataset</th><th>Weight</th></tr></thead>
          <tbody>${r.datasets.map((d, i) => `
            <tr><td class="num">${i + 1}</td>
              <td><a href="#/dataset/${d.id}">${esc(d.title)}</a>
                <br><span class="sub">${esc([d.first_author, d.year].filter(Boolean).join(' '))}
                ${d.geo_id ? ' · ' + esc(d.geo_id) : ''} · ${d.num_conds} conditions</span></td>
              <td class="bar"><span>${d.weight.toFixed(1)}%</span>
                <i style="width:${(100 * d.weight / maxWeight).toFixed(1)}%"></i></td></tr>`).join('')}
          </tbody></table></div>
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
      body: JSON.stringify({ gene_ids: geneIds, aspect, max_results: 15 }),
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
  box.innerHTML = `<div class="table-wrap"><table>
    <thead><tr><th>GO term</th><th>In list</th><th>Genome</th><th>Corrected p</th><th>FDR</th><th>Genes</th></tr></thead>
    <tbody>${sig.map((t) => `<tr>
      <td><a href="https://www.ebi.ac.uk/QuickGO/term/${encodeURIComponent(t.term_id)}" target="_blank"
        rel="noopener">${esc(t.name)}</a><br><span class="sub nowrap">${esc(t.term_id)}</span></td>
      <td class="num">${t.list_count}/${t.list_size}</td>
      <td class="num">${t.background_count}/${t.background_size}</td>
      <td class="num">${p(t.p_bonferroni)}</td>
      <td class="num">${p(t.fdr)}</td>
      <td class="desc"><div title="${esc(t.genes.join(', '))}">${t.genes.map(esc).join(', ')}</div></td>
    </tr>`).join('')}</tbody></table></div>
    <p class="sub">Hypergeometric test, Bonferroni-corrected p and Benjamini-Hochberg FDR. ${source}</p>`;
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

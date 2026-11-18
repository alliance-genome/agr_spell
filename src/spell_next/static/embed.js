'use strict';

// Compact single-gene result for iframes: /embed?gene=ACT1&genes=10&sets=5

(async () => {
  const box = document.getElementById('embed');
  const params = new URLSearchParams(location.search);
  const gene = (params.get('gene') || '').trim();
  const nGenes = Math.min(50, Math.max(1, Number(params.get('genes')) || 10));
  const nSets = Math.min(20, Math.max(0, Number(params.get('sets') ?? 5)));
  const esc = (s) => String(s ?? '').replace(/[&<>"']/g,
    (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const label = (g) => g.standard || g.systematic;
  const full = (hash) => `${location.origin}/${hash}`;
  const fail = (msg) => { box.innerHTML = `<p class="notice">${esc(msg)}</p>`; };

  if (!gene) return fail('No gene given (use ?gene=NAME).');
  try {
    const lookup = await (await fetch(`/api/genes/lookup?q=${encodeURIComponent(gene)}`)).json();
    if (!lookup.genes || lookup.genes.length !== 1) return fail(`Gene not found or ambiguous: ${gene}`);
    const res = await fetch(`/api/genes/${lookup.genes[0].id}/coexpression?limit=${nGenes}` +
      `&datasets=${Math.max(1, nSets)}&via=embed`);
    if (!res.ok) return fail('Result not available.');
    const r = await res.json();
    const q = r.query[0];
    box.innerHTML = `
      <h1>Genes co-expressed with ${esc(label(q))}</h1>
      <ol class="embed-list">${r.genes.map((g, i) => `<li><span class="rank">${i + 1}</span>
        <a class="name" href="${full(`#/search/${encodeURIComponent(g.systematic)}`)}" target="_blank"
          rel="noopener" title="${esc(g.description)}">${esc(label(g))}</a>
        <span class="score">${g.score.toFixed(2)}</span></li>`).join('')}</ol>
      ${nSets ? `<h2>Most informative datasets</h2>
      <ol class="embed-list">${r.datasets.slice(0, nSets).map((d, i) => `<li><span class="rank">${i + 1}</span>
        <a class="name" href="${full(`#/dataset/${d.id}`)}" target="_blank" rel="noopener"
          title="${esc(d.title)}">${esc(d.title)}</a>
        <span class="score">${d.weight.toFixed(1)}%</span></li>`).join('')}</ol>` : ''}
      <a class="more" href="${full(`#/search/${encodeURIComponent(q.systematic)}`)}" target="_blank"
        rel="noopener">Full results in SPELL ↗</a>`;
  } catch (err) {
    fail(err.message);
  }
})();

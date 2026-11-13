"""GO term enrichment for SPELL result gene lists (replaces SPELL's GO::TermFinder service).

Annotations are propagated to all ancestors over is_a and part_of. Each term with at
least two list genes is tested with the hypergeometric distribution against all genes
annotated in the same aspect, with Bonferroni and Benjamini-Hochberg corrections.

Parent terms often contain exactly the same list genes as a more significant term, so by
default a term whose list genes match a higher-ranked term is folded into that term's
``also`` list instead of being reported separately. Corrections still count every test.
"""
import gzip
import math
from collections import defaultdict
from dataclasses import dataclass, field

ASPECTS = {'P': 'biological_process', 'F': 'molecular_function', 'C': 'cellular_component'}
NAMESPACE_ASPECT = {v: k for k, v in ASPECTS.items()}


def read_obo(path):
    """Return {term id: (name, aspect, [parent ids])} for non-obsolete terms."""
    terms, cur = {}, None

    def flush():
        if cur and cur.get('id') and not cur.get('obsolete') and cur.get('aspect'):
            terms[cur['id']] = (cur.get('name', cur['id']), cur['aspect'], cur['parents'])

    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('['):
                flush()
                cur = {'parents': []} if line == '[Term]' else None
            elif cur is not None and ': ' in line:
                key, val = line.split(': ', 1)
                if key == 'id':
                    cur['id'] = val
                elif key == 'name':
                    cur['name'] = val
                elif key == 'namespace':
                    cur['aspect'] = NAMESPACE_ASPECT.get(val)
                elif key == 'is_obsolete' and val == 'true':
                    cur['obsolete'] = True
                elif key == 'is_a':
                    cur['parents'].append(val.split(' ', 1)[0])
                elif key == 'relationship' and val.startswith('part_of '):
                    cur['parents'].append(val.split(' ')[1])
    flush()
    return terms


def ancestors(terms):
    """Return {term: set of the term and all its ancestors}."""
    memo = {}

    def walk(t):
        if t not in memo:
            memo[t] = {t}  # placeholder guards against cycles
            acc = {t}
            for p in terms.get(t, (None, None, []))[2]:
                if p in terms:
                    acc |= walk(p)
            memo[t] = acc
        return memo[t]

    for t in terms:
        walk(t)
    return memo


def read_gaf(path, sgdid_to_gene, terms):
    """Return {gene id: set of directly annotated term ids} (NOT annotations skipped)."""
    out = defaultdict(set)
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('!'):
                continue
            f = line.split('\t')
            if len(f) < 9 or 'NOT' in f[3]:
                continue
            gene = sgdid_to_gene.get(f[1])
            if gene is not None and f[4] in terms:
                out[gene].add(f[4])
    return out


def propagate(direct, anc):
    return {g: set().union(*(anc[t] for t in ts)) for g, ts in direct.items()}


def log_comb(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k, M, K, N):
    """P(X >= k) for X ~ Hypergeometric(population M, K successes, N draws)."""
    hi = min(K, N)
    if k > hi:
        return 0.0
    base = log_comb(M, N)
    logs = [log_comb(K, i) + log_comb(M - K, N - i) - base for i in range(k, hi + 1)]
    top = max(logs)
    return min(1.0, math.exp(top) * sum(math.exp(x - top) for x in logs))


@dataclass
class Enrichment:
    term_id: str
    name: str
    aspect: str
    list_count: int
    list_size: int
    background_count: int
    background_size: int
    p_value: float
    p_bonferroni: float
    fdr: float
    genes: list
    also: list = field(default_factory=list)   # [(term id, name)] with the same list genes


class GOIndex:
    """In-memory GO annotations for the compendium's genes."""

    def __init__(self, term_info, gene_terms):
        self.term_info = term_info                      # term -> (name, aspect)
        self.gene_terms = gene_terms                    # gene -> set of terms (propagated)
        self.term_genes = defaultdict(set)
        for g, ts in gene_terms.items():
            for t in ts:
                self.term_genes[t].add(g)
        self.background = {a: {g for g, ts in gene_terms.items()
                               if any(self.term_info[t][1] == a for t in ts)} for a in ASPECTS}

    def enrich(self, genes, aspect='P', min_list_count=2, max_results=50, collapse=True):
        bg = self.background[aspect]
        listed = [g for g in dict.fromkeys(genes) if g in bg]
        M, N = len(bg), len(listed)
        if N == 0:
            return []
        counts = defaultdict(list)
        for g in listed:
            for t in self.gene_terms[g]:
                if self.term_info[t][1] == aspect:
                    counts[t].append(g)
        tested = []
        for t, hits in counts.items():
            if len(hits) < min_list_count:
                continue
            K = len(self.term_genes[t])
            if K == M:              # root terms carry no information
                continue
            tested.append((hypergeom_sf(len(hits), M, K, N), t, hits, K))
        n_tests = len(tested)
        tested.sort()
        out, prev_fdr = [], 1.0
        fdrs = [0.0] * n_tests
        for i in range(n_tests - 1, -1, -1):        # Benjamini-Hochberg, monotone from the top
            prev_fdr = min(prev_fdr, tested[i][0] * n_tests / (i + 1))
            fdrs[i] = prev_fdr
        by_genes = {}
        for (p, t, hits, K), fdr in zip(tested, fdrs):
            name, asp = self.term_info[t]
            key = frozenset(hits)
            if collapse and key in by_genes:
                by_genes[key].also.append((t, name))
                continue
            if len(out) >= max_results:
                continue        # keep scanning so shown rows still collect their duplicates
            row = Enrichment(t, name, asp, len(hits), N, K, M, p, min(1.0, p * n_tests), fdr, hits)
            by_genes[key] = row
            out.append(row)
        return out

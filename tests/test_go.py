import math
from itertools import combinations

import pytest

from spell_next import go

OBO = """format-version: 1.2
data-version: releases/test

[Term]
id: GO:0000001
name: root process
namespace: biological_process

[Term]
id: GO:0000002
name: child process
namespace: biological_process
is_a: GO:0000001 ! root process

[Term]
id: GO:0000003
name: part of child
namespace: biological_process
relationship: part_of GO:0000002 ! child process

[Term]
id: GO:0000004
name: obsolete thing
namespace: biological_process
is_obsolete: true

[Typedef]
id: part_of
name: part of
"""


def test_hypergeom_sf_matches_brute_force():
    M, K, N = 12, 5, 4
    population = [1] * K + [0] * (M - K)
    draws = list(combinations(range(M), N))
    for k in range(N + 1):
        brute = sum(sum(population[i] for i in d) >= k for d in draws) / len(draws)
        assert go.hypergeom_sf(k, M, K, N) == pytest.approx(brute)
    assert go.hypergeom_sf(5, M, K, N) == 0.0


def test_read_obo_and_propagation(tmp_path):
    path = tmp_path / 'go.obo'
    path.write_text(OBO)
    terms = go.read_obo(str(path))
    assert set(terms) == {'GO:0000001', 'GO:0000002', 'GO:0000003'}
    anc = go.ancestors(terms)
    assert anc['GO:0000003'] == {'GO:0000001', 'GO:0000002', 'GO:0000003'}
    assert go.propagate({7: {'GO:0000003'}}, anc)[7] == anc['GO:0000003']


def _index():
    info = {'GO:R': ('root', 'P'), 'GO:A': ('alpha', 'P'), 'GO:B': ('beta', 'P'), 'GO:F': ('fun', 'F')}
    gene_terms = {g: {'GO:R'} for g in range(100)}
    for g in range(10):
        gene_terms[g] |= {'GO:A'}
    for g in range(10, 40):
        gene_terms[g] |= {'GO:B'}
    gene_terms[0] |= {'GO:F'}
    return go.GOIndex(info, gene_terms)


def test_enrich_finds_the_planted_term_and_skips_root():
    res = _index().enrich(list(range(8)) + [50, 60], 'P')
    assert [r.term_id for r in res] == ['GO:A']
    top = res[0]
    assert (top.list_count, top.list_size, top.background_count, top.background_size) == (8, 10, 10, 100)
    assert top.p_value == pytest.approx(go.hypergeom_sf(8, 100, 10, 10))
    assert top.p_bonferroni == pytest.approx(top.p_value)    # one term tested


def test_enrich_corrections_with_two_terms():
    res = _index().enrich([0, 1, 2, 10, 11], 'P')
    assert {r.term_id for r in res} == {'GO:A', 'GO:B'}
    for r in res:
        assert r.p_bonferroni == pytest.approx(min(1.0, 2 * r.p_value))
        assert r.p_value <= r.fdr <= r.p_bonferroni + 1e-12
    assert math.isclose(res[-1].fdr, res[-1].p_value)        # largest p keeps its own value


def test_enrich_ignores_genes_outside_the_aspect_background():
    assert _index().enrich([500, 501], 'P') == []


def test_terms_with_identical_list_genes_are_collapsed():
    info = {'GO:R': ('root', 'P'), 'GO:A': ('alpha', 'P'), 'GO:PA': ('parent of alpha', 'P'),
            'GO:B': ('beta', 'P')}
    gene_terms = {g: {'GO:R'} for g in range(100)}
    for g in range(10):
        gene_terms[g] |= {'GO:A', 'GO:PA'}
    for g in range(10, 20):
        gene_terms[g] |= {'GO:PA'}          # the parent is bigger in the genome ...
    for g in range(20, 30):
        gene_terms[g] |= {'GO:B'}
    index = go.GOIndex(info, gene_terms)
    genes = list(range(6)) + [20, 21, 22]  # ... but holds the same list genes as GO:A
    collapsed = index.enrich(genes, 'P')
    assert [r.term_id for r in collapsed] == ['GO:A', 'GO:B']
    assert collapsed[0].also == [('GO:PA', 'parent of alpha')]
    full = index.enrich(genes, 'P', collapse=False)
    assert [r.term_id for r in full] == ['GO:A', 'GO:PA', 'GO:B']
    assert collapsed[0].p_bonferroni == full[0].p_bonferroni        # corrections count all 3 tests

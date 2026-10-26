from pathlib import Path

from openacme.icd10 import Icd10Graph
from openacme.acme import get_acme_graph, get_acme_tables, make_valid_range

RESOURCES = Path(__file__).parent / 'resources'
ACME_FILE = str(RESOURCES / 'acme_test.html')
ICD10_FILE = str(RESOURCES / 'icd10_test.xml.zip')


def test_make_valid_range():
    g = Icd10Graph()
    # These are all real ranges that are problematic from the ACME tables
    assert make_valid_range('B59', 'B64', g) == ('B60', 'B64')
    assert make_valid_range('I84.0', 'I84.9', g) is None
    assert make_valid_range('D76.0', 'D86.9', g) == ('D76.1', 'D86.9')
    assert make_valid_range('M31.2', 'M31.3', g) == 'M31.3'
    assert make_valid_range('U01.0', 'U02', g) is None
    assert make_valid_range('A91', 'A92.4', g) == ('A92', 'A92.4')
    assert make_valid_range('B58.0', 'B59', g) == ('B58.0', 'B58.9')


def test_acme_node_kinds():
    g = get_acme_graph(acme_file=ACME_FILE, icd10_file=ICD10_FILE)
    assert len(g.nodes) > 0
    for n, data in g.nodes(data=True):
        if isinstance(n, tuple):
            assert data['kind'] == 'range'
        else:
            assert data['kind'] == 'category'


def test_acme_range_expansion_nodes_present():
    icd10 = Icd10Graph(icd10_file=ICD10_FILE)
    g = get_acme_graph(acme_file=ACME_FILE, icd10_file=ICD10_FILE)
    ranges = [n for n, d in g.nodes(data=True) if d['kind'] == 'range']
    assert ranges
    for start, end in ranges:
        expanded = icd10.expand_icd10_range(start, end)
        assert expanded
        for code in expanded:
            assert code in g
            assert g.nodes[code]['kind'] == 'category'


def test_acme_tables():
    tables = get_acme_tables(acme_file=ACME_FILE, icd10_file=ICD10_FILE)
    g = tables['graph']
    assert g.edges[('A00.0', 'A01.1'), 'A00.1']['maybe']
    # One row is marked 'M' and the other is not
    assert not g.edges[('A01.0', 'A01.1'), 'A00.0']['maybe']
    assert tables['table_e'] == [
        {'target': 'A00.1', 'source': 'A00.0', 'rule_type': 'DS',
         'result': None, 'maybe': False, 'hash': False, 'reason': None},
        {'target': 'A00.1', 'source': ('A01.0', 'A01.1'), 'rule_type': 'DSC',
         'result': 'A00.9', 'maybe': True, 'hash': True,
         'reason': 'Suba must be acute'},
    ]
    assert tables['code_lists'] == {
        'A': ['A00.0', 'A01'],
        'B': ['A00.9'],
        'C1': ['A01.1'],
        'C2': ['A01.0'],
        'C3': ['A00', 'A00.0', 'A00.1', 'A00.9', 'A01', 'A01.0', 'A01.1'],
        'H': ['A00.1', 'A01', 'A01.0', 'A01.1'],
    }

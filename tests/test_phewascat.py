from ontoma import OnToma

def test_answer():
    t = OnToma()
    assert t.efo_lookup('asthma') == 'EFO:0000270'
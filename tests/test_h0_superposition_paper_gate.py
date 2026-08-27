import numpy as np


def test_delta_mu_roundtrip():
    from tools.h0_superposition_paper_gate import delta_mu, h_from_delta_mu
    h1=71.92804067892604; h2=73.0
    dm=delta_mu(h2,h1)
    assert abs(h_from_delta_mu(h1,dm)-h2)<1e-10


def test_subset_rows_collapses_repeated_photometry():
    from tools.h0_superposition_paper_gate import subset_rows
    keys=np.array(['N5584_2007af','N5584_2007af','N1365_2012fr',''])
    got=subset_rows(keys,['N5584_2007af'])
    assert got.tolist()==[0,1]


def test_gate_ratios_are_dimensionless():
    from tools.h0_superposition_paper_gate import residual_fraction
    assert abs(residual_fraction(0.6,1.2)-0.5)<1e-12

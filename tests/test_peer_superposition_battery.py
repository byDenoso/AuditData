import math


def test_equiv_dmu_is_exact_inverse_h0_ratio():
    from tools.peer_superposition_battery import equiv_dmu
    h1=71.92804067892604
    h2=73.0
    d=equiv_dmu(h2,h1)
    assert abs(d - (-5*math.log10(h2/h1))) < 1e-12


def test_combined_local_sigma_includes_peer_and_hole_envelope():
    from tools.peer_superposition_battery import combined_local_sigma
    s=combined_local_sigma(
        peer_h0=70.391,
        peer_sigma=0.801,
        local_median=71.92804067892604,
        local_lo=71.71203503729637,
        local_hi=72.17037753571651,
    )
    expected_peer=0.801*(71.92804067892604/70.391)
    expected_hole=((72.17037753571651-71.71203503729637)/2)/math.sqrt(3)
    assert abs(s-math.sqrt(expected_peer**2+expected_hole**2)) < 1e-12


def test_rank_signed_influence_prefers_largest_downward_shift():
    from tools.peer_superposition_battery import rank_signed_influence
    rows=[
        {'id':'a','delta_H0':0.4},
        {'id':'b','delta_H0':-0.2},
        {'id':'c','delta_H0':-0.7},
    ]
    r=rank_signed_influence(rows,direction='down')
    assert [x['id'] for x in r]==['c','b','a']
    r=rank_signed_influence(rows,direction='up')
    assert [x['id'] for x in r]==['a','b','c']

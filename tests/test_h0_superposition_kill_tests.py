import numpy as np


def test_canonical_cid_strips_shoes_suffix_but_not_core_name():
    from tools.h0_superposition_kill_tests import canonical_cid
    assert canonical_cid('1992bl_50') == '1992bl'
    assert canonical_cid('PSNJ0252467_150') == 'psnj0252467'
    assert canonical_cid('ASASSN-16hh_18') == 'asassn16hh'
    assert canonical_cid('2005df_ANU') == '2005dfanu'


def test_redshift_y_shift_is_zero_for_identical_redshifts_and_sign_is_physical():
    from tools.h0_superposition_kill_tests import redshift_y_shift
    assert abs(redshift_y_shift(0.03, 0.03)) < 1e-14
    # Replacing a larger flow-corrected z by a smaller raw z makes y less negative,
    # which lowers inferred H0 for y ~ MB - 5logH0.
    assert redshift_y_shift(0.03, 0.025) > 0


def test_deleted_estimator_weight_reproduces_cached_fit():
    from tools.h0_superposition_kill_tests import build_precision_cache, deleted_estimator_weight, fit_delete_cached
    y=np.array([1.1,1.9,3.2,3.9,5.1])
    L=np.column_stack([np.ones(5), np.arange(5.)])
    C=np.array([[1,.1,0,0,0],[.1,1,.1,0,0],[0,.1,1,.1,0],[0,0,.1,1,.1],[0,0,0,.1,1.]])
    params=['a','5logH0']
    cache=build_precision_cache(y,L,C)
    D=np.array([1,3])
    q,cov,p=fit_delete_cached(y,L,params,cache,D,())
    w=deleted_estimator_weight(L,params,cache,D,(),target='5logH0')
    assert abs(w@y - q[p.index('5logH0')]) < 1e-10


def test_mw_anchor_rows_are_explicit_geometric_constraints_only():
    from tools.h0_superposition_kill_tests import anchor_rows
    src=np.array(['MHW1_HST','MHW1_Gaia','MW','LMC_HST','mu_LMC','mu_N4258'])
    assert anchor_rows(src,'MW').tolist() == [0,1]
    assert anchor_rows(src,'LMC').tolist() == [4]
    assert anchor_rows(src,'N4258').tolist() == [5]

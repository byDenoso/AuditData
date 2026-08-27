import numpy as np


def test_gls_recovers_linear_solution():
    from tools.shoes_host_jackknife import gls_fit
    y=np.array([1.,2.,3.,4.])
    L=np.array([[1.,0.],[1.,1.],[1.,2.],[1.,3.]])
    C=np.eye(4)
    q,cov=gls_fit(y,L,C)
    assert np.allclose(q,[1.,1.],atol=1e-10)
    assert cov.shape==(2,2)


def test_drop_host_removes_exact_and_calibrator_rows_and_parameter():
    from tools.shoes_host_jackknife import drop_host
    sources=np.array(['N1015','N1015','N1015_2009ig_51','OTHER'])
    params=['mu_N1015','MB','5logH0']
    y=np.arange(4.,dtype=float)
    L=np.arange(12.,dtype=float).reshape(4,3)
    C=np.eye(4)
    yy,LL,CC,pp=drop_host(y,L,C,sources,params,'N1015')
    assert list(pp)==['MB','5logH0']
    assert yy.tolist()==[3.]
    assert LL.shape==(1,2)
    assert CC.shape==(1,1)


def test_h0_from_fivelogh0():
    from tools.shoes_host_jackknife import h0_from_fivelogh0
    assert abs(h0_from_fivelogh0(5*np.log10(73.04))-73.04)<1e-10


def test_anchor_constraint_rows_remove_only_geometric_prior():
    from tools.shoes_host_jackknife import anchor_constraint_rows
    s=np.array(['MHW1_HST','MHW1_Gaia','LMC_HST','LMC_GRND','SMC','mu_LMC','N4258','mu_N4258'])
    assert anchor_constraint_rows(s,'LMC').tolist()==[5]
    assert anchor_constraint_rows(s,'N4258').tolist()==[7]
    assert anchor_constraint_rows(s,'MW').tolist()==[0,1]


def test_superposition_canonical_cid_and_redshift_shift():
    from tools.h0_superposition_kill_tests import canonical_cid, redshift_y_shift
    assert canonical_cid('1992bl_50')=='1992bl'
    assert canonical_cid('PSNJ0252467_150')=='psnj0252467'
    assert canonical_cid('ASASSN-16hh_18')=='asassn16hh'
    assert canonical_cid('2005df_ANU')=='2005dfanu'
    assert abs(redshift_y_shift(.03,.03))<1e-14
    assert redshift_y_shift(.03,.025)>0


def test_deleted_estimator_weight_reproduces_cached_fit():
    from tools.shoes_host_jackknife import build_precision_cache, fit_delete_cached
    from tools.h0_superposition_kill_tests import deleted_estimator_weight
    y=np.array([1.1,1.9,3.2,3.9,5.1])
    L=np.column_stack([np.ones(5),np.arange(5.)])
    C=np.array([[1,.1,0,0,0],[.1,1,.1,0,0],[0,.1,1,.1,0],[0,0,.1,1,.1],[0,0,0,.1,1.]])
    params=['a','5logH0']
    cache=build_precision_cache(y,L,C)
    D=np.array([1,3])
    q,cov,p=fit_delete_cached(y,L,params,cache,D,())
    w=deleted_estimator_weight(L,params,cache,D,(),target='5logH0')
    assert abs(w@y-q[p.index('5logH0')])<1e-10

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

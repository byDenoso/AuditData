import numpy as np


def test_tail_ranks_identify_low_density_and_high_flow():
    from tools.local_joint_control import tail_ranks
    x = np.array([10.,20.,30.,40.])
    low = tail_ranks(x, lower=True)
    high = tail_ranks(x, lower=False)
    assert low[0] < low[-1]
    assert high[-1] < high[0]
    assert np.all((low > 0) & (low <= 1))
    assert np.all((high > 0) & (high <= 1))


def test_joint_score_prefers_low_density_high_outflow():
    from tools.local_joint_control import joint_score
    pd = np.array([0.01, 0.5, 0.5])
    pv = np.array([0.02, 0.5, 0.01])
    s = joint_score(pd, pv)
    assert s[0] > s[1]
    assert s[2] > s[1]


def test_rotation_maps_unit_vectors_to_unit_vectors():
    from tools.local_joint_control import random_rotation
    rng = np.random.default_rng(9)
    R = random_rotation(rng)
    assert np.allclose(R.T @ R, np.eye(3), atol=1e-12)
    assert np.isclose(np.linalg.det(R), 1.0, atol=1e-12)

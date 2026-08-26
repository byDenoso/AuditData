import numpy as np


def test_fibonacci_directions_are_unit_vectors():
    from tools.local_universe_nearby_scan import fibonacci_sphere
    u = fibonacci_sphere(64)
    assert u.shape == (64, 3)
    assert np.allclose(np.linalg.norm(u, axis=1), 1.0, atol=1e-12)


def test_shuffle_positions_preserves_radii():
    from tools.local_universe_nearby_scan import shuffle_positions
    rng = np.random.default_rng(7)
    pos = rng.normal(size=(200, 3))
    pos *= rng.uniform(10, 200, size=(200, 1)) / np.linalg.norm(pos, axis=1)[:, None]
    shuffled = shuffle_positions(pos, np.random.default_rng(8))
    assert np.allclose(np.sort(np.linalg.norm(pos, axis=1)), np.sort(np.linalg.norm(shuffled, axis=1)))


def test_scan_finds_inserted_void_near_true_center():
    from tools.local_universe_nearby_scan import scan_counts_against_shuffles
    rng = np.random.default_rng(11)
    n = 6000
    u = rng.normal(size=(n, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
    r = 180 * rng.random(n) ** (1/3)
    pos = u * r[:, None]
    true_center = np.array([60.0, 0.0, 0.0])
    true_R = 35.0
    pos = pos[np.linalg.norm(pos - true_center, axis=1) > true_R]
    centers = np.array([
        true_center,
        [0.0, 60.0, 0.0],
        [0.0, 0.0, 60.0],
        [-60.0, 0.0, 0.0],
    ])
    out = scan_counts_against_shuffles(pos, centers, radius=true_R, nnull=24, seed=13)
    best = int(np.argmin(out['z']))
    assert best == 0
    assert out['ratio'][0] < 0.6


def test_outflow_regression_recovers_positive_slope():
    from tools.local_universe_nearby_scan import outflow_regression
    rng = np.random.default_rng(17)
    center = np.array([45.0, -10.0, 5.0])
    n = 1200
    u = rng.normal(size=(n, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
    rr = rng.uniform(25, 110, n)
    pos = center + u * rr[:, None]
    nhat = pos / np.linalg.norm(pos, axis=1)[:, None]
    x = np.sum((pos-center) * nhat, axis=1)
    vpec = 3.2 * x + rng.normal(0, 120, n)
    fit = outflow_regression(pos, vpec, center)
    assert fit['slope'] > 2.0
    assert fit['n'] == n

import numpy as np


def test_required_calibration_shift_matches_exact_distance_modulus_relation():
    from tools.peer_h0_paper_battery import required_delta_mu
    # PEER+Local-Hole central observer prediction from the frozen prior battery.
    got = required_delta_mu(71.92804067892604, 73.0)
    assert np.isclose(got, -5*np.log10(73.0/71.92804067892604), atol=1e-12)


def test_profiled_common_offset_is_zero_when_prediction_equals_data():
    from tools.peer_h0_paper_battery import profile_common_offset
    obs = np.array([72.0, 72.0, 72.0])
    sig = np.array([1.0, 1.5, 2.0])
    fit = profile_common_offset(obs, sig, predicted_h0=72.0)
    assert abs(fit['delta_mu']) < 1e-12
    assert abs(fit['delta_h0']) < 1e-12


def test_leaveout_gate_requires_survival_without_dominant_calibrator():
    from tools.peer_h0_paper_battery import leaveout_gate
    # A paper claim must fail if its full residual disappears after dropping one calibrator.
    result = leaveout_gate(required_residual=1.0, shifts=np.array([1.05, 0.2, 0.1]), max_fraction=0.75)
    assert result['passes'] is False
    assert result['max_fraction_of_residual'] > 1.0


def test_effects_are_not_declared_additive_without_covariance():
    from tools.peer_h0_paper_battery import combine_effects
    out = combine_effects(np.array([0.60, 0.53]), covariance=None)
    assert out['identifiable'] is False
    assert np.isnan(out['combined_shift'])


def test_mock_null_returns_calibrated_tail_probability():
    from tools.peer_h0_paper_battery import empirical_pvalue
    null = np.array([0.1,0.2,0.3,0.4,0.5])
    # plus-one rule: (1 + 2 values >= .4) / (1+5) = 0.5
    assert np.isclose(empirical_pvalue(0.4, null, tail='greater'), 0.5)

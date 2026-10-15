"""Correlated-series statistics, against synthetic data with known answers.

An AR(1) process has an analytic correlation time, so these are verification
tests: the estimator is checked against mathematics, not against another
estimator.
"""

from __future__ import annotations

import numpy as np
import pytest

from saxsnmr.stats import (
    assess_convergence,
    autocorrelation_time,
    block_average,
    equilibration_start,
)


def ar1(n: int, phi: float, seed: int = 0) -> np.ndarray:
    """First-order autoregressive series with correlation coefficient phi.

    Its integrated autocorrelation time is tau = 0.5 * (1 + phi) / (1 - phi).
    """
    rng = np.random.default_rng(seed)
    noise = rng.normal(size=n)
    x = np.empty(n)
    x[0] = noise[0]
    for i in range(1, n):
        x[i] = phi * x[i - 1] + noise[i]
    return x


def tau_exact(phi: float) -> float:
    return 0.5 * (1.0 + phi) / (1.0 - phi)


def test_white_noise_has_no_correlation():
    """Independent samples give tau = 0.5, the no-correlation floor."""
    rng = np.random.default_rng(0)
    tau = autocorrelation_time(rng.normal(size=20000))
    np.testing.assert_allclose(tau, 0.5, atol=0.35)


@pytest.mark.parametrize("phi", [0.5, 0.8, 0.9])
def test_autocorrelation_time_matches_the_ar1_value(phi):
    """Checked against the closed form, not against another implementation."""
    tau = autocorrelation_time(ar1(60000, phi, seed=1))
    np.testing.assert_allclose(tau, tau_exact(phi), rtol=0.25)


def test_correlated_data_inflates_the_error_bar():
    """The whole point: pretending independence understates the error.

    For phi = 0.9 the correlation time is 9.5 samples, so the honest error
    is roughly sqrt(19) ~ 4.4 times the naive one. Reporting the naive value
    is how an unconverged run looks precise.
    """
    stats = block_average(ar1(20000, 0.9, seed=2))
    assert stats.sem > stats.naive_sem
    np.testing.assert_allclose(stats.inflation, np.sqrt(2 * tau_exact(0.9)), rtol=0.3)
    assert stats.n_effective < stats.n_samples / 10


def test_independent_data_is_not_inflated():
    stats = block_average(np.random.default_rng(3).normal(size=5000))
    np.testing.assert_allclose(stats.inflation, 1.0, atol=0.3)


def test_block_average_refuses_a_single_sample():
    """A one-frame ensemble has no spread; returning 0 would invent one."""
    with pytest.raises(ValueError, match="at least 2 samples"):
        block_average(np.array([1.0]))


def test_equilibration_is_detected_and_discarded():
    """A decaying transient followed by a stationary tail."""
    n_burn, n_prod = 400, 4000
    burn = 12.0 * np.exp(-np.arange(n_burn) / 60.0)
    series = np.concatenate([burn, np.zeros(n_prod)])
    series = series + np.random.default_rng(4).normal(scale=0.2, size=series.size)
    start = equilibration_start(series)
    assert start > 100, f"burn-in not detected (start={start})"
    assert start < n_burn + 400, f"discarded far too much (start={start})"


def test_a_drifting_series_is_reported_as_unconverged():
    """A mean taken from a drifting series describes no particular state."""
    n = 2000
    drifting = np.linspace(0, 6, n) + np.random.default_rng(5).normal(scale=0.1, size=n)
    report = assess_convergence(drifting, "Rg")
    assert not report.converged
    assert "drifts" in report.reason


def test_a_stationary_series_is_reported_as_converged():
    stationary = np.random.default_rng(6).normal(size=4000)
    report = assess_convergence(stationary, "Rg")
    assert report.converged, report.reason


def test_a_short_correlated_series_is_unresolved_not_converged():
    """Ten frames of a slowly varying quantity cannot resolve a mean.

    This is the case the MD notebooks are actually in: a 5 ps run sampled
    every 0.5 ps gives ten frames of a correlated observable.
    """
    report = assess_convergence(ar1(10, 0.9, seed=7), "Rg")
    assert not report.converged
    assert "too short" in report.reason

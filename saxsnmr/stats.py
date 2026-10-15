"""Error bars for correlated time series.

Successive MD frames are not independent. Treating them as if they were makes
the standard error too small by the square root of the correlation time,
which is how a trajectory that has not converged can report a tight, confident
and wrong number.

Everything here is a verification tool rather than physics: the functions have
answers that can be checked against synthetic data with a known correlation
structure, and the tests do exactly that.

References
----------
Flyvbjerg & Petersen, J. Chem. Phys. 91, 461 (1989), block averaging.
    https://doi.org/10.1063/1.457480
Chodera, J. Chem. Theory Comput. 12, 1799 (2016), equilibration detection.
    https://doi.org/10.1021/acs.jctc.5b00784
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "BlockAverage",
    "ConvergenceReport",
    "assess_convergence",
    "autocorrelation_time",
    "block_average",
    "equilibration_start",
    "statistical_inefficiency",
]


def autocorrelation_time(x: np.ndarray) -> float:
    """Integrated autocorrelation time, in samples.

    Summed with the initial-positive-sequence truncation: the estimator of
    rho(t) becomes pure noise once the signal has decayed, and summing that
    noise to the end of the series adds variance without adding information.
    A value near 0.5 means successive samples are effectively independent.
    """
    x = np.asarray(x, dtype=float).ravel()
    n = x.size
    if n < 4:
        return float("nan")
    centred = x - x.mean()
    var = float(centred @ centred) / n
    if var <= 0:
        return 0.5  # a constant series has no correlation
    tau = 0.5
    for lag in range(1, n // 2):
        rho = float(centred[:-lag] @ centred[lag:]) / ((n - lag) * var)
        if rho <= 0.0:
            break
        tau += rho
    return float(tau)


def statistical_inefficiency(x: np.ndarray) -> float:
    """g = 2 tau: how many samples one independent observation is worth.

    The effective sample size is len(x) / g, and the standard error of the
    mean is sqrt(g) times the naive one.
    """
    tau = autocorrelation_time(x)
    return float("nan") if not np.isfinite(tau) else max(2.0 * tau, 1.0)


@dataclass
class BlockAverage:
    """Mean of a correlated series with an honest uncertainty."""

    mean: float
    sem: float  # correlation-corrected standard error
    naive_sem: float  # what pretending independence would give
    n_samples: int
    n_effective: float
    statistical_inefficiency: float

    @property
    def inflation(self) -> float:
        """How much the naive error bar understates the real one."""
        return self.sem / self.naive_sem if self.naive_sem > 0 else float("nan")

    def __repr__(self) -> str:
        return (
            f"BlockAverage({self.mean:.3f} +/- {self.sem:.3f}, "
            f"n_eff={self.n_effective:.1f} of {self.n_samples}, "
            f"g={self.statistical_inefficiency:.1f})"
        )


def block_average(x: np.ndarray) -> BlockAverage:
    """Mean and correlation-corrected standard error of a time series.

    Raises
    ------
    ValueError
        If fewer than two samples are given. A single frame has no spread,
        and returning zero or NaN would let a one-frame "ensemble" be
        reported with an error bar.
    """
    x = np.asarray(x, dtype=float).ravel()
    if x.size < 2:
        raise ValueError(f"need at least 2 samples to estimate an error, got {x.size}")
    g = statistical_inefficiency(x)
    naive = float(x.std(ddof=1) / np.sqrt(x.size))
    n_eff = x.size / g if g > 0 else float(x.size)
    return BlockAverage(
        mean=float(x.mean()),
        sem=float(x.std(ddof=1) / np.sqrt(max(n_eff, 1.0))),
        naive_sem=naive,
        n_samples=int(x.size),
        n_effective=float(n_eff),
        statistical_inefficiency=float(g),
    )


def equilibration_start(x: np.ndarray) -> int:
    """Index after which the series looks equilibrated.

    Chooses the truncation that maximises the number of *effective* samples
    in what remains: discarding burn-in removes drift and shortens the
    correlation time, and the gain stops once only equilibrated data is left.
    Returns 0 when discarding nothing is best.
    """
    x = np.asarray(x, dtype=float).ravel()
    n = x.size
    if n < 8:
        return 0
    best_t, best_eff = 0, -np.inf
    for t in range(0, n - 4, max(1, n // 50)):
        tail = x[t:]
        g = statistical_inefficiency(tail)
        eff = tail.size / g if np.isfinite(g) and g > 0 else tail.size
        if eff > best_eff:
            best_t, best_eff = t, eff
    return int(best_t)


@dataclass
class ConvergenceReport:
    """Whether a trajectory supports the number taken from it."""

    quantity: str
    equilibration_index: int
    production: BlockAverage
    drift_per_sample: float
    drift_over_run: float
    converged: bool
    reason: str

    def __repr__(self) -> str:
        verdict = "converged" if self.converged else "NOT converged"
        return (
            f"ConvergenceReport({self.quantity}: {self.production.mean:.3f} "
            f"+/- {self.production.sem:.3f}, {verdict} -- {self.reason})"
        )


def assess_convergence(
    x: np.ndarray, quantity: str = "value", min_effective: float = 10.0
) -> ConvergenceReport:
    """Judge whether a time series has equilibrated and been sampled enough.

    Two separate failures are checked. A series can drift steadily, in which
    case its mean describes no particular state; or it can be stationary but
    so short that the mean is unresolved. Both are reported, because the fix
    differs: longer equilibration against a longer production run.

    ``min_effective`` is the number of independent samples below which the
    mean is called unresolved. Ten is a weak bar, chosen so that failing it
    is unambiguous rather than marginal.
    """
    x = np.asarray(x, dtype=float).ravel()
    start = equilibration_start(x)
    tail = x[start:]
    if tail.size < 2:
        raise ValueError("no samples remain after discarding equilibration")
    stats = block_average(tail)

    # Least-squares slope over the production part, tested against the
    # uncertainty *of the slope* -- not against the error on the mean.
    #
    # Those differ by a large constant factor. For an uncorrelated series the
    # fitted drift across the window has a standard deviation of about
    # sqrt(12) times the standard error of the mean, purely from fitting a
    # line to noise. Comparing against the mean's error therefore flags
    # stationary data as drifting most of the time, which is useless: a
    # convergence check that always cries wolf gets ignored.
    idx = np.arange(tail.size, dtype=float)
    slope = intercept = 0.0
    se_slope = float("inf")
    if tail.size > 3:
        slope, intercept = (float(v) for v in np.polyfit(idx, tail, 1))
        resid = tail - (slope * idx + intercept)
        sxx = float(((idx - idx.mean()) ** 2).sum())
        if sxx > 0:
            s = float(np.sqrt((resid @ resid) / max(tail.size - 2, 1)))
            # Inflate by sqrt(g): correlated residuals carry less information
            # about the slope than their count suggests.
            se_slope = (
                s / np.sqrt(sxx) * np.sqrt(max(stats.statistical_inefficiency, 1.0))
            )
    drift = slope * tail.size

    reasons = []
    if stats.n_effective < min_effective:
        reasons.append(
            f"only {stats.n_effective:.1f} effective samples "
            f"(g={stats.statistical_inefficiency:.1f}); the run is too short"
        )
    if np.isfinite(se_slope) and se_slope > 0 and abs(slope) > 2.0 * se_slope:
        reasons.append(
            f"drifts by {drift:+.3f} over the production window "
            f"({abs(slope) / se_slope:.1f} sigma on the slope)"
        )
    return ConvergenceReport(
        quantity=quantity,
        equilibration_index=start,
        production=stats,
        drift_per_sample=slope,
        drift_over_run=drift,
        converged=not reasons,
        reason="; ".join(reasons) if reasons else "stationary and sampled",
    )

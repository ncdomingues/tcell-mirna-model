# -*- coding: utf-8 -*-
"""Small statistics helpers -- no scipy dependency."""
import math


def wilson_ci(k, n, z=1.96):
    """Wilson score interval for a binomial proportion. Returns (pct, lo, hi)
    as percentages. More reliable than the naive normal approximation when
    the proportion is near 0 or 100 (which several nodes here are)."""
    if n == 0:
        return 0.0, 0.0, 0.0
    p = k / n
    denom = 1 + z**2 / n
    centre = p + z**2 / (2 * n)
    adj = z * math.sqrt((p * (1 - p) + z**2 / (4 * n)) / n)
    lo = (centre - adj) / denom
    hi = (centre + adj) / denom
    return round(100 * p, 1), round(100 * max(0, lo), 1), round(100 * min(1, hi), 1)

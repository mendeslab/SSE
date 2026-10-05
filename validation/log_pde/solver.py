"""Experimental one-dimensional log transport with analytic exterior values."""
import math
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import logsumexp


# Each mixture row is (log integrated weight, mean, variance). Evolution of each
# Gaussian is exact for backward drift v*d/dx and variance rate q on the real line.
def log_reference(x, components, t=0.0, drift=0.0, diffusion=0.2):
    weights, means, variances = np.asarray(components).T
    variances = variances + diffusion * t
    means = means - drift * t
    terms = (weights[:, None] - .5 * np.log(2 * np.pi * variances[:, None])
             - (np.asarray(x)[None, :] - means[:, None]) ** 2 / (2 * variances[:, None]))
    return logsumexp(terms, axis=0)


# Advance mixture parameters analytically, retaining integrated component weights.
def advance(components, t, drift, diffusion):
    result = np.array(components, dtype=float, copy=True)
    result[:, 1] -= drift * t
    result[:, 2] += diffusion * t
    return result


# Multiply every pair of normalized Gaussian components and include the overlap
# integral in its weight. Keep all pairs: merging/truncation would change the reference.
def merge(left, right, birth=0.15):
    result = []
    for w1, m1, s1 in left:
        for w2, m2, s2 in right:
            variance = s1 * s2 / (s1 + s2)
            mean = (m1 * s2 + m2 * s1) / (s1 + s2)
            overlap = -.5 * (math.log(2 * math.pi * (s1 + s2)) + (m1 - m2) ** 2 / (s1 + s2))
            result.append((math.log(birth) + w1 + w2 + overlap, mean, variance))
    return np.array(result)


# Centered weights reproduce derivatives of polynomials through the stated order.
# The first derivative uses antisymmetric pairs, the second symmetric pairs; these
# explicit formulas avoid a separate coefficient-generation algorithm in the solver.
def derivatives(values, dx, order):
    first = {2: [1 / 2], 4: [2 / 3, -1 / 12], 6: [3 / 4, -3 / 20, 1 / 60]}[order]
    second = {2: [-2, 1], 4: [-5 / 2, 4 / 3, -1 / 12],
              6: [-49 / 18, 3 / 2, -3 / 20, 1 / 90]}[order]
    radius = order // 2
    size = len(values) - 2 * radius
    center = values[radius:radius + size]
    slope = np.zeros(size)
    # Difference from center avoids cancellation of a large spatially uniform offset.
    curvature = np.zeros(size)
    for j in range(1, radius + 1):
        below = values[radius - j:radius - j + size]
        above = values[radius + j:radius + j + size]
        slope += first[j - 1] * (above - below)
        curvature += second[j] * ((above - center) + (below - center))
    return slope / dx, curvature / dx ** 2


# Solve h_t = v*h_x + q/2*(h_xx+h_x^2) at all stored points, using exact
# time-dependent reference values outside the grid. A fixed log offset is restored
# on output. It has zero spatial/time derivatives, so does not alter the equation.
def propagate(x, initial, components, duration, drift, diffusion, order, atol, max_evaluations=20000):
    dx = x[1] - x[0]
    radius = order // 2
    exterior = np.r_[x[0] - dx * np.arange(radius, 0, -1), x[-1] + dx * np.arange(1, radius + 1)]
    offset = float(np.max(initial))
    evaluations = 0
    started = time.perf_counter()

    # Exterior values come from the same full analytic problem at each RK stage.
    # Reject nonfinite evolution or exhausted work explicitly; neither is an accuracy result.
    def rhs(t, h):
        nonlocal evaluations
        evaluations += 1
        if evaluations > max_evaluations:
            raise RuntimeError("derivative-evaluation budget exceeded")
        outside = log_reference(exterior, components, t, drift, diffusion) - offset
        padded = np.r_[outside[:radius], h, outside[radius:]]
        with np.errstate(over="ignore", invalid="ignore"):
            slope, curvature = derivatives(padded, dx, order)
            result = drift * slope + diffusion / 2 * (curvature + slope ** 2)
        if not np.all(np.isfinite(result)):
            raise RuntimeError("nonfinite log derivative")
        return result

    try:
        solution = solve_ivp(rhs, (0, duration), initial - offset, method="RK45", atol=atol, rtol=1e-12)
        status = "ok" if solution.success else solution.message
        result = solution.y[:, -1] + offset if solution.success else None
        times = solution.t
    except RuntimeError as error:
        status, result, times = str(error), None, np.array([])
    return result, dict(status=status, evaluations=evaluations, seconds=time.perf_counter() - started,
                        times=times)

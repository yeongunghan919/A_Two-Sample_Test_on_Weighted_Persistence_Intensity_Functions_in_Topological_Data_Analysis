import numpy as np
import scipy.spatial
import math
import pandas as pd
from itertools import combinations
from persim import wasserstein
from persim.landscapes import (
    PersLandscapeApprox,
    average_approx,
    snap_pl,
    plot_landscape,
    plot_landscape_simple
)
import random
from ripser import ripser
from itertools import product
from scipy.stats import norm

"""
================================================================================
Bonferroni test for equality of intensity functions over a collection of bandwidths
H0:p=q vs H1:p\neq q
================================================================================
"""
# ============================================================
#  Bandwidth-aggregated Bonferroni test
# ============================================================

def estimate_bandwidth_reference(
    sample_diagrams,
    max_point_pairs=2_000_000
):
    """Estimate (r_b, r_d) from cross-diagram coordinate gaps.

    Take coordinate-wise medians of absolute differences between points in
    *different* diagrams. Each pair of nonempty diagrams receives equal total
    weight, even when their point counts differ. With two points per diagram,
    this equals the ordinary median of all cross-diagram point-pair gaps used
    in the practical-bandwidth simulation. When there are more than
    max_point_pairs cross-diagram point pairs, sample an equal number of point
    pairs within each diagram pair, with a fixed seed for reproducibility.
    """
    if int(max_point_pairs) != max_point_pairs or max_point_pairs < 1:
        raise ValueError("max_point_pairs must be a positive integer.")
    max_point_pairs = int(max_point_pairs)

    diagrams = []
    for diagram in sample_diagrams:
        diagram = np.asarray(diagram, dtype=float)
        if diagram.size == 0:
            diagram = np.empty((0, 2), dtype=float)
        if diagram.ndim != 2 or diagram.shape[1] != 2:
            raise ValueError("Each diagram must have shape (points, 2).")
        if not np.all(np.isfinite(diagram)):
            raise ValueError("Point coordinates must be finite.")
        if len(diagram):
            diagrams.append(diagram)

    if len(diagrams) < 2:
        raise ValueError("At least two nonempty diagrams are required.")

    counts = [len(diagram) for diagram in diagrams]
    total_pairs = sum(counts[i] * sum(counts[i + 1:])
                      for i in range(len(counts) - 1))
    diagram_pairs = len(diagrams) * (len(diagrams) - 1) // 2
    if diagram_pairs > max_point_pairs:
        raise ValueError(
            "max_point_pairs must allow at least one point pair per "
            "pair of nonempty diagrams."
        )

    balanced = len(set(counts)) == 1
    subsample = total_pairs > max_point_pairs
    if balanced and not subsample:
        points = np.vstack(diagrams)
        owners = np.repeat(np.arange(len(diagrams)), counts[0])
        cross_diagram = owners[:, None] < owners[None, :]
        differences = np.abs(points[:, None, :] - points[None, :, :])
        reference = np.median(differences[cross_diagram], axis=0)
    else:
        gaps = []
        pair_weights = []
        rng = np.random.default_rng(0)
        per_diagram_pair = max_point_pairs // diagram_pairs
        for i, first in enumerate(diagrams):
            for second in diagrams[i + 1:]:
                count = len(first) * len(second)
                if subsample and count > per_diagram_pair:
                    first_ids = rng.integers(len(first), size=per_diagram_pair)
                    second_ids = rng.integers(len(second), size=per_diagram_pair)
                    difference = np.abs(first[first_ids] - second[second_ids])
                else:
                    difference = np.abs(first[:, None, :] - second[None, :, :])
                    difference = difference.reshape(-1, 2)
                gaps.append(difference)
                if not balanced:
                    # Total weight one per pair of diagrams, independently
                    # of diagram size or the number of sampled point pairs.
                    pair_weights.append(np.full(len(difference), 1.0 / len(difference)))

        gaps = np.concatenate(gaps, axis=0)
        if balanced:
            reference = np.median(gaps, axis=0)
        else:
            weights = np.concatenate(pair_weights)
            reference = np.empty(2, dtype=float)
            for coordinate in range(2):
                order = np.argsort(gaps[:, coordinate], kind="stable")
                cumulative = np.cumsum(weights[order])
                midpoint = 0.5 * cumulative[-1]
                position = int(np.searchsorted(cumulative, midpoint))
                value = gaps[order[position], coordinate]
                # Agree with np.median when exactly half the weight lies below.
                if (position + 1 < len(order)
                        and np.isclose(cumulative[position], midpoint,
                                       rtol=1e-12, atol=0.0)):
                    value = 0.5 * (value + gaps[order[position + 1], coordinate])
                reference[coordinate] = value

    if np.any(reference <= 0.0):
        raise ValueError("Both coordinate medians must be positive.")
    return reference


def practical_bandwidth_collection(n, m, reference_scales):
    """Return the nine (2**j r_b/sqrt(N), 2**k r_d/sqrt(N)) pairs."""
    reference = np.asarray(reference_scales, dtype=float)
    if reference.shape != (2,) or not np.all(np.isfinite(reference)):
        raise ValueError("reference_scales must be two finite numbers.")
    if np.any(reference <= 0.0):
        raise ValueError("Both reference_scales must be positive.")
    if n + m <= 0:
        raise ValueError("n+m must be positive.")

    base = reference / np.sqrt(n + m)
    return [
        (float(2.0**j * base[0]), float(2.0**k * base[1]))
        for j, k in product((-1, 0, 1), repeat=2)
    ]


def bonferroni_test(
    X,
    Y,
    alpha=0.05,
    func_weight=None,
    reference_fraction=0.20,
    reference_rng=None,
    max_block_entries=1_000_000,
):
    """
    Apply the studentized Bonferroni test with nine practical bandwidths.

    Compute coordinate-wise reference scales from the current X and Y sample.
    The bandwidth collection therefore changes with each test sample.

    Parameters
    ----------
    X, Y : list of np.ndarray
        Samples of persistence diagrams.
    alpha : float
        Significance level.
    func_weight : callable or None
        A function taking one point (birth, death) and returning
        its weight. If None, the constant weight is used.
    max_block_entries : int
        Maximum number of point pairs processed in one memory block.
        Reduce this value if memory is especially limited.

    Returns
    -------
    result : dict
        Test decision and diagnostic information.
    """
    X = list(X)
    Y = list(Y)

    n = len(X)
    m = len(Y)

    if n < 4 or m < 4:
        raise ValueError(
            "The variance estimators require n >= 4 and m >= 4."
        )

    if not (0.0 < alpha < 1.0):
        raise ValueError("alpha must lie in (0, 1).")

    # --------------------------------------------------------
    # Data-dependent coordinate-wise scales and nine practical bandwidths
    # --------------------------------------------------------
    if not 0 < reference_fraction <= 1:
        raise ValueError("reference_fraction must lie in (0, 1].")

    rng = (
        np.random.default_rng(0)
        if reference_rng is None
        else reference_rng
    )

    n_ref = max(1, int(np.ceil(reference_fraction * n)))
    m_ref = max(1, int(np.ceil(reference_fraction * m)))

    reference_X = [
        X[i] for i in rng.choice(n, size=n_ref, replace=False)
    ]
    reference_Y = [
        Y[j] for j in rng.choice(m, size=m_ref, replace=False)
    ]

    reference_scales = estimate_bandwidth_reference(
        reference_X + reference_Y
    )
    bandwidths = practical_bandwidth_collection(
        n, m, reference_scales
    )
    # reference_scales = estimate_bandwidth_reference(X + Y)
    # bandwidths = practical_bandwidth_collection(n, m, reference_scales)

    number_of_bandwidths = len(bandwidths)

    critical_value = float(
        norm.isf(alpha / number_of_bandwidths)
    )

    # --------------------------------------------------------
    # Statistics over all bandwidths
    # --------------------------------------------------------
    geometry = prepare_geometry(
        X + Y,
        func_weight=func_weight
    )

    # Compute all Gram matrices in memory-bounded blocks.  In particular,
    # this avoids materialising several total_points x total_points arrays.
    grams = diagram_grams(
        geometry,
        bandwidths,
        max_block_entries=max_block_entries
    )

    statistics = []
    u_statistics = []
    variance_estimates = []

    for gram in grams:
        statistic, u_statistic, sigma_squared = T_statistic(
            gram,
            n,
            m
        )

        statistics.append(statistic)
        u_statistics.append(u_statistic)
        variance_estimates.append(sigma_squared)

    statistics = np.asarray(statistics, dtype=float)
    u_statistics = np.asarray(u_statistics, dtype=float)
    variance_estimates = np.asarray(
        variance_estimates,
        dtype=float
    )

    valid = (
        np.isfinite(statistics)
        & np.isfinite(variance_estimates)
        & (variance_estimates > 0.0)
    )

    if np.any(valid):
        best_index = int(np.argmax(statistics))
        maximum_statistic = float(statistics[best_index])
        best_bandwidth = bandwidths[best_index]

        adjusted_p_value = min(
            1.0,
            number_of_bandwidths
            * norm.sf(maximum_statistic)
        )
    else:
        maximum_statistic = -np.inf
        best_bandwidth = None
        adjusted_p_value = 1.0

    reject = maximum_statistic > critical_value

    return {
        "reject": bool(reject),
        "adjusted_p_value": float(adjusted_p_value),
        # "maximum_statistic": maximum_statistic,
        # "critical_value": critical_value,
        # "best_bandwidth": best_bandwidth,
        # "number_of_bandwidths": number_of_bandwidths,
        # "number_of_invalid_variances": int(np.sum(~valid)),
         "bandwidths": bandwidths
        # "statistics": statistics,
        # "u_statistics": u_statistics,
        # "variance_estimates": variance_estimates,
    }
# ============================================================
# 1. Precompute diagram geometry
# ============================================================

def prepare_geometry(diagrams, func_weight=None):
    """
    Precompute quantities shared across all bandwidths.

    Unlike the previous implementation, this function does not construct
    total_points x total_points distance and index matrices.  Its memory use
    is linear in the total number of persistence points.

    Parameters
    ----------
    diagrams : list of np.ndarray
        Each diagram has shape (number_of_points, 2).
    func_weight : callable or None
        A function taking one point (birth, death) and returning
        its weight. If None, constant weight is used.
    """
    if func_weight is None:
        func_weight = function_weight("constant")

    cleaned_diagrams = []

    for diagram in diagrams:
        diagram = np.asarray(diagram, dtype=float)

        if diagram.size == 0:
            diagram = np.empty((0, 2), dtype=float)

        if diagram.ndim != 2 or diagram.shape[1] != 2:
            raise ValueError(
                "Each persistence diagram must have shape "
                "(number_of_points, 2)."
            )

        cleaned_diagrams.append(diagram)

    number_of_diagrams = len(cleaned_diagrams)

    counts = np.array(
        [len(diagram) for diagram in cleaned_diagrams],
        dtype=int
    )

    offsets = np.concatenate((
        np.array([0], dtype=np.int64),
        np.cumsum(counts, dtype=np.int64)
    ))

    total_points = int(offsets[-1])

    if total_points == 0:
        points = np.empty((0, 2), dtype=float)
    else:
        points = np.vstack([
            diagram
            for diagram in cleaned_diagrams
            if len(diagram) > 0
        ])

    owners = np.repeat(
        np.arange(number_of_diagrams, dtype=np.int64),
        counts
    )

    # Apply func_weight to every persistence point
    weights = np.fromiter(
        (func_weight(point) for point in points),
        dtype=float,
        count=len(points)
    )

    if not np.all(np.isfinite(weights)):
        raise ValueError(
            "func_weight returned a non-finite value."
        )

    return (
        points,
        weights,
        owners,
        offsets,
        number_of_diagrams
    )

# ============================================================
# 2. Diagram-kernel Gram matrix
# ============================================================

def diagram_grams(
    geometry,
    bandwidths,
    max_block_entries=1_000_000
):
    """
    Construct diagram-level Gram matrices for all bandwidths using
    memory-bounded blocks.

        k_lambda(x,y)
          = exp[-0.5 * {
                (x1-y1)^2/lambda_1^2
                + (x2-y2)^2/lambda_2^2
            }]
            / (2*pi*lambda_1*lambda_2)

    The result follows the arithmetic and accumulation order of the original
    full-matrix implementation, but processes consecutive rows in blocks.
    The largest point-pair array contains at most approximately
    ``max_block_entries`` floating-point values.

    Parameters
    ----------
    geometry : tuple
        Output of ``prepare_geometry``.
    bandwidths : sequence of pairs
        Each pair is (lambda_1, lambda_2).
    max_block_entries : int
        Maximum number of point pairs in a processing block.  The default
        keeps each float64 block near 8 MB.

    Returns
    -------
    grams : np.ndarray
        Array of shape (number_of_bandwidths, number_of_diagrams,
        number_of_diagrams).
    """
    bandwidths = np.asarray(bandwidths, dtype=float)

    if bandwidths.ndim != 2 or bandwidths.shape[1] != 2:
        raise ValueError(
            "bandwidths must have shape (number_of_bandwidths, 2)."
        )

    if len(bandwidths) == 0:
        raise ValueError("At least one bandwidth is required.")

    if not np.all(np.isfinite(bandwidths)) or np.any(bandwidths <= 0.0):
        raise ValueError("Bandwidths must be finite and positive.")

    if int(max_block_entries) != max_block_entries or max_block_entries < 1:
        raise ValueError("max_block_entries must be a positive integer.")

    max_block_entries = int(max_block_entries)

    (
        points,
        weights,
        owners,
        offsets,
        number_of_diagrams
    ) = geometry

    number_of_bandwidths = len(bandwidths)
    grams = np.zeros(
        (
            number_of_bandwidths,
            number_of_diagrams,
            number_of_diagrams
        ),
        dtype=float
    )

    if len(points) == 0:
        return grams

    total_points = len(points)
    rows_per_block = max(
        1,
        max_block_entries // total_points
    )

    for diagram_index in range(number_of_diagrams):
        diagram_start = int(offsets[diagram_index])
        diagram_stop = int(offsets[diagram_index + 1])

        if diagram_start == diagram_stop:
            continue

        one_block_for_diagram = (
            diagram_stop - diagram_start <= rows_per_block
        )

        for row_start in range(
            diagram_start,
            diagram_stop,
            rows_per_block
        ):
            row_stop = min(
                row_start + rows_per_block,
                diagram_stop
            )
            row_points = points[row_start:row_stop]

            squared_birth_differences = (
                row_points[:, None, 0]
                - points[None, :, 0]
            )**2

            squared_death_differences = (
                row_points[:, None, 1]
                - points[None, :, 1]
            )**2

            row_weights = weights[row_start:row_stop]
            weight_products = (
                row_weights[:, None] * weights[None, :]
            )

            # These are the diagram labels of the flattened block in C
            # (row-major) order, matching diagram_pair_ids.ravel() in the
            # original implementation.
            block_owners = np.tile(
                owners,
                row_stop - row_start
            )

            for bandwidth_index, (
                lambda_1,
                lambda_2
            ) in enumerate(bandwidths):
                # Preserve the original order of floating-point operations:
                # division by lambda**2, addition, multiplication by -0.5,
                # exponential, normalization, and finally weight products.
                weighted_point_kernel = (
                    squared_birth_differences
                    / lambda_1**2
                )
                weighted_point_kernel += (
                    squared_death_differences
                    / lambda_2**2
                )
                weighted_point_kernel *= -0.5
                np.exp(
                    weighted_point_kernel,
                    out=weighted_point_kernel
                )
                weighted_point_kernel /= (
                    2.0 * np.pi * lambda_1 * lambda_2
                )
                weighted_point_kernel *= weight_products

                flattened_kernel = weighted_point_kernel.ravel()

                if one_block_for_diagram:
                    # A single bincount has exactly the same input order as
                    # the original full-matrix implementation.
                    grams[
                        bandwidth_index,
                        diagram_index,
                        :
                    ] = np.bincount(
                        block_owners,
                        weights=flattened_kernel,
                        minlength=number_of_diagrams
                    )
                else:
                    # np.add.at accumulates repeated indices sequentially.
                    # Consequently, splitting a large diagram across blocks
                    # does not change the original summation order.
                    np.add.at(
                        grams[
                            bandwidth_index,
                            diagram_index,
                            :
                        ],
                        block_owners,
                        flattened_kernel
                    )

    return grams


def diagram_gram(
    geometry,
    lambda_1,
    lambda_2,
    max_block_entries=1_000_000
):
    """Construct one memory-bounded diagram-level Gram matrix."""
    return diagram_grams(
        geometry,
        [(lambda_1, lambda_2)],
        max_block_entries=max_block_entries
    )[0]


# ============================================================
# 3. Variance-component estimators
# ============================================================

def a_hat(gram):
    """
    Unbiased estimator of tr(Sigma^2).

    This is the O(r^2) implementation of the fourth-order
    U-statistic estimator.
    """
    gram = np.asarray(gram, dtype=float)
    r = gram.shape[0]

    if gram.shape != (r, r):
        raise ValueError("gram must be a square matrix.")

    if r < 4:
        raise ValueError("At least four diagrams are required.")

    off_diagonal = gram.copy()
    np.fill_diagonal(off_diagonal, 0.0)

    row_sums = off_diagonal.sum(axis=1)
    total_sum = row_sums.sum()

    numerator = (
        total_sum**2
        - 2.0 * (r - 1) * np.sum(row_sums**2)
        + (r - 1) * (r - 2)
        * np.sum(off_diagonal**2)
    )

    denominator = r * (r - 1) * (r - 2) * (r - 3)

    return float(numerator / denominator)


def a_pq_hat(gram_xy):
    """
    Unbiased estimator of tr(Sigma_P Sigma_Q).
    """
    gram_xy = np.asarray(gram_xy, dtype=float)
    n, m = gram_xy.shape

    if n < 2 or m < 2:
        raise ValueError(
            "gram_xy must have at least two rows and columns."
        )

    row_sums = gram_xy.sum(axis=1)
    column_sums = gram_xy.sum(axis=0)
    total_sum = gram_xy.sum()

    numerator = (
        total_sum**2
        - n * np.sum(row_sums**2)
        - m * np.sum(column_sums**2)
        + n * m * np.sum(gram_xy**2)
    )

    denominator = n * (n - 1) * m * (m - 1)

    return float(numerator / denominator)


# ============================================================
# 4. Studentized statistic
# ============================================================

def T_statistic(gram, n, m):
    """
    Compute the studentized statistic.

    Returns
    -------
    statistic : float
        Studentized test statistic.
    u_statistic : float
        Unstudentized U-statistic.
    sigma_squared : float
        Estimated null variance.
    """
    gram = np.asarray(gram, dtype=float)

    if gram.shape != (n + m, n + m):
        raise ValueError(
            "gram must have shape (n+m, n+m)."
        )

    gram_xx = gram[:n, :n]
    gram_yy = gram[n:, n:]
    gram_xy = gram[:n, n:]

    sum_xx_off_diagonal = (
        gram_xx.sum() - np.trace(gram_xx)
    )

    sum_yy_off_diagonal = (
        gram_yy.sum() - np.trace(gram_yy)
    )

    u_statistic = (
        sum_xx_off_diagonal / (n * (n - 1))
        + sum_yy_off_diagonal / (m * (m - 1))
        - 2.0 * gram_xy.sum() / (n * m)
    )

    a_p = a_hat(gram_xx)
    a_q = a_hat(gram_yy)
    a_pq = a_pq_hat(gram_xy)

    sigma_squared = (
        2.0 * a_p / (n * (n - 1))
        + 2.0 * a_q / (m * (m - 1))
        + 4.0 * a_pq / (n * m)
    )

    if (
        not np.isfinite(sigma_squared)
        or sigma_squared <= 0.0
    ):
        statistic = -np.inf
    else:
        statistic = u_statistic / np.sqrt(sigma_squared)

    return (
        float(statistic),
        float(u_statistic),
        float(sigma_squared)
    )




"""
================================================================================
Aggregated permutation test for intensity functions over a collection of bandwidths
H0:P=Q vs H1:p\neq q
This code is adapted from the code of Schrab et al., 2023.
================================================================================
"""

def Aggtest(
    X,
    Y,
    alpha=0.05,
    kernel="gaussian", # Only the Gaussian kernel is defined in this code. If you want a different kernel, you must input it yourself.
    number_bandwidths=10,
    optimal_bandwidths=True,
    weight_function=None,
    Rff_approx=False,
    B=1000,
    merging_function=np.min,
    seed=42,
    return_dictionary=False,
):
    """
    Two sample Aggtest using persistence diagrams.
    
    Given data from one distribution and data from another distribution,
    return 0 if the test fails to reject the null 
    (i.e. data comes from the same distribution), 
    or return 1 if the test rejects the null 
    (i.e. their intensity functions are different).
    
    Parameters
    ----------
    X: list of persistence diagrams
    Y: list of persistence diagrams
    alpha: scalar
        The value of alpha must be between 0 and 1.
    kernel: str
        The value of kernel must be "gaussian".
        Only the Gaussian kernel is defined in this code. If you want a different kernel, you must input it yourself.
    number_bandwidths: int
        The number of bandwidths per kernel to include in the collection.
    weight_function: function
        function to compute the weight of each point in persistence diagram.
    B: int
        The number of Monte Carlo samples of permutations to generate.
    seed: int 
        Random seed used for the randomness of the permutations.
    return_dictionary: bool
        If true, a dictionary is returned containing for step-by-step results.
    Returns
    -------
    output : int
        0 if the Aggtest fails to reject the null
            (i.e. data comes from the same distribution)
        1 if the Aggtest rejects the null
            (i.e. intensity functions are different)
    dictionary: dict
        Returned only if return_dictionary is True.
        Dictionary containing the overall output of the Aggtest.

    """    
    # Assertions
    n = len(X)
    m = len(Y)
    assert 0 < alpha  and alpha < 1
    assert kernel in (
        "gaussian",
    )
    assert number_bandwidths > 1 and type(number_bandwidths) == int
    assert B>0 and type(B)==int

    if isinstance(weight_function, (list, tuple)):
        weight_functions = list(weight_function)
    else:
        weight_functions = [weight_function]

    assert len(weight_functions) > 0
    assert all(callable(w) for w in weight_functions)

    if not optimal_bandwidths:
        # Coordinate-wise median-based bandwidth collections
        def compute_bandwidths(distances, number_bandwidths):
            if np.min(distances) < 10 ** (-1):
                d = np.sort(distances)
                lambda_min = np.maximum(
                    d[int(np.floor(len(d) * 0.05))],
                    10 ** (-1)
                )
            else:
                lambda_min = np.min(distances)
            lambda_min = lambda_min / 2
            lambda_max = np.maximum(
                np.max(distances),
                3 * 10 ** (-1)
            )
            lambda_max = lambda_max * 2
            power = (lambda_max / lambda_min) ** (
                1 / (number_bandwidths - 1)
            )
            bandwidths = np.array([
                power ** i * lambda_min
                for i in range(number_bandwidths)
            ])
            return bandwidths
        max_samples = 100
        all_dists_1 = []
        all_dists_2 = []
        for diagramX in X[:max_samples]:
            for diagramY in Y[:max_samples]:
                # Skip empty diagrams
                if len(diagramX) == 0 or len(diagramY) == 0:
                    continue
                # Coordinate 1 pairwise distances
                d1 = np.abs(
                    diagramX[:, 0][:, None]
                    - diagramY[:, 0][None, :]
                ).reshape(-1)
                # Coordinate 2 pairwise distances
                d2 = np.abs(
                    diagramX[:, 1][:, None]
                    - diagramY[:, 1][None, :]
                ).reshape(-1)

                all_dists_1.append(d1)
                all_dists_2.append(d2)
        distances_1 = np.concatenate(all_dists_1)
        distances_2 = np.concatenate(all_dists_2)
        bandwidths_1 = compute_bandwidths(
            distances_1,
            number_bandwidths
        )
        bandwidths_2 = compute_bandwidths(
            distances_2,
            number_bandwidths
        )
        # Cartesian product:
        # (lambda_1, lambda_2)
        bandwidths = np.array(
            list(product(bandwidths_1, bandwidths_2)),
            dtype=float
        )
        number_bandwidths = len(bandwidths)


    else: #optimal bandwidth collection which ensures the minimax optimality up to an iterated logarithm factor.
        T = math.ceil(math.log2((n + m) / math.log(math.log(n + m))))

        bandwidth_1d = np.array([2**(-i) for i in range(1, T + 1)])

        bandwidths = np.array(list(product(bandwidth_1d, repeat=2)))

        number_bandwidths = len(bandwidths)

    # Setup permutations   #Efficient permutation method (Schrab et al., 2023, p. 51)
    rs = np.random.RandomState(seed)
    idx = rs.rand(B, n+m).argsort(axis=1)  # (B, n+m): rows of permuted indices
    #11
    v11 = np.concatenate((np.ones(n), -np.ones(m)))  # (n+m, )
    V11i = np.tile(v11, (B, 1))  # (B, n+m)
    V11 = np.take_along_axis(V11i, idx, axis=1)  # (B, n+m): permute the entries of the rows
    V11 = np.vstack([V11, v11])    # (B+1)th entry is the original test statistic (no permutation)
    V11 = V11.transpose()  # (n+m, B+1)
    #10
    v10 = np.concatenate((np.ones(n), np.zeros(m)))
    V10i = np.tile(v10, (B, 1))
    V10 = np.take_along_axis(V10i, idx, axis=1)
    V10 = np.vstack([V10, v10])
    V10 = V10.transpose()
    #01
    v01 = np.concatenate((np.zeros(n), -np.ones(m)))
    V01i = np.tile(v01, (B, 1))
    V01 = np.take_along_axis(V01i, idx, axis=1)
    V01= np.vstack([V01, v01])
    V01 = V01.transpose()

        
    # Step 1: compute the statistical matrix M
    parameter_grid = list(product(bandwidths, weight_functions))

    N = len(parameter_grid)
    M = np.zeros((N, B + 1))
    list_pd = X + Y

    for i, (bandwidth, weight) in enumerate(parameter_grid):
        K = Mat_gram(
            list_pd,
            bandwidth,
            weight,
            kernel,
            seed=seed,
            Rff_approx=Rff_approx,
            num_rff=10**4,
        )

        # Set diagonal elements to zero
        np.fill_diagonal(K, 0)

        # Compute permuted test statistics
        M[i] = (
            np.sum(V10 * (K @ V10), 0)
            * (m - n + 1) / (m * n * (n - 1))
            + np.sum(V01 * (K @ V01), 0)
            * (n - m + 1) / (m * n * (m - 1))
            + np.sum(V11 * (K @ V11), 0) / (m * n)
        )

    M = M.transpose()

    # Step 2: compute P-value Matrix
    def computing_rank(vector):
        vector = np.asarray(vector)
        L = len(vector)
        return np.sum(vector[:, None] <= vector[None, :], axis=1) / L
    Pvar_M=np.column_stack([computing_rank(M[:, j]) for j in range(M.shape[1])]) #P-value matrix

    # Step 3: Aggregated Statistics
    Agg_Stats=np.apply_along_axis(merging_function,axis=1,arr=Pvar_M)
    # Step 4: output test result
    pval= np.sum(Agg_Stats<=Agg_Stats[-1])/len(Agg_Stats) #Agg_Stats[-1] corresponds to the identity permutation.

    
    # create rejection dictionary 
    reject_dictionary = {}

    reject_dictionary["Bandwidth"] = np.array(
        [bandwidth for bandwidth, _ in parameter_grid],
         dtype=float)
    reject_dictionary["Weight function"] = [
             weight for _, weight in parameter_grid
    ]
    reject_dictionary["Stat_Matrix"] = M
    reject_dictionary["P-value_Matrix"] = Pvar_M
    reject_dictionary["AggStats"] = Agg_Stats
    reject_dictionary["p-value"] = pval
    # Aggregated test rejects if pval is less than or equal to alpha
    reject_dictionary["Aggtest reject"] = True if pval <= alpha else False

    if return_dictionary:
        return int(reject_dictionary["Aggtest reject"]), reject_dictionary
    else:
        return int(reject_dictionary["Aggtest reject"])



def gaussian_kernel_func(x, y, bandwidth):
    bandwidth = np.asarray(bandwidth)
    return np.exp(-0.5 * np.sum(((x - y) / bandwidth)**2))

def Linear(diagram_1, diagram_2, vec_weight_1, vec_weight_2,kernel,bandwidth):
    """
    diagram_1, _2: Input diagrams, np.arraies of shape (n,2)
    vec_weight_1 , _2: Lists of length n consisting of weight values of points in diagram_1, _2, respectively.
    kernel : A string representing the kernel function
    bandwidth : Flaot, a bandwidth of kernel function
    
    return : Float, a linear kernel value between diagram_1 and diagram_2
    """
    if kernel=='gaussian':
        kernel_func = gaussian_kernel_func
    
    return_value = 0.0
    num_point_1 = diagram_1.shape[0] 
    num_point_2 = diagram_2.shape[0]
    for i in range(num_point_1):
        for j in range(num_point_2):
            return_value += (vec_weight_1[i] * vec_weight_2[j]
                  * kernel_func(diagram_1[i, :], diagram_2[j, :],bandwidth))
    return return_value

def function_weight(name_weight, arc_c=1.0, arc_p=1.0, lin_el=1.0,poly_order=2):
    if name_weight == "arctan":
        def func_weight(bd):
            return np.maximum(np.arctan(math.pow((bd[1] - bd[0]), arc_p) / arc_c), 0.0)
    elif name_weight == "linear":
        def func_weight(bd):
            return (bd[1]-bd[0])
    elif name_weight == "Poly":
        def func_weight(bd):
            return (bd[1]-bd[0])**poly_order
    else:  # unweighted
        def func_weight(bd):
            return 1.0
    return func_weight

def vector_weight(function_weight,diagram):
        num_point = diagram.shape[0]
        vec = np.empty(num_point)
        for k in range(num_point):
            vec[k] = function_weight(diagram[k, :])
        return vec

def list_vector_weight(list_pd,function_weight):
    list_weight = []
    for i in range(len(list_pd)):
        list_weight.append(vector_weight(function_weight=function_weight,diagram=list_pd[i]))
    return list_weight

def Mat_gram(list_pd,bandwidth,weight_function,kernel,seed,Rff_approx=False,num_rff=10**4):
    """
    list_pd: List of length (n+m) which contains all input persistence diagrams
    kernel : String representing the kernel function
    bandwidth : a list representing the bandwidth of kernel function
    weight_function : Function, the weight_function for computing weight of each point in each diagram.
    Rff_approx : Boolean, 1 -> Random Fourier Feature approximation for computing the Linear-kernel value between two diagrams.
    num_rff : The number of samples for the Rff approximation
    
    return : An (n+m) * (n+m) matrix whose element is the Linear kernel value between the corresponding two diagrams.
    """
    num_pd = len(list_pd)
    weight_values = list_vector_weight(list_pd,function_weight= weight_function)
    random.seed(seed)
    np.random.seed(seed)
    if Rff_approx and kernel=='gaussian':
        mat_rff = np.empty((num_pd, num_rff))
        Z= np.random.multivariate_normal(
                [0.0, 0.0], [[bandwidth[0] ** (-2.0), 0.0], [0.0, bandwidth[1] ** (-2.0)]], num_rff)
        b= np.random.uniform(0.0,2.0*math.pi,num_rff)
        for k in range(num_pd):
            mat_rff[k,:]= np.dot(weight_values[k],np.sqrt(2)*np.cos(np.inner(list_pd[k],Z)+b))
        mat =  np.inner(mat_rff, mat_rff) / num_rff
    else:
        mat = np.empty((num_pd, num_pd))
        for i in range(num_pd):
            for j in range(i + 1):
                mat[i, j] = Linear(
                    list_pd[i], list_pd[j],
                    weight_values[i], weight_values[j],kernel,bandwidth)
                mat[j, i] = mat[i, j]
                
    return mat



"""
================================================================================
Persistence Diagram permutation test

This implementation follows Robinson and Turner, 2017.
================================================================================
"""

def compute_distance_matrix(pd_list):
    """
    Computes symmetric matrix of all pairwise Wasserstein distances
    """
    n = len(pd_list)
    dist_mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = wasserstein(pd_list[i], pd_list[j])
            dist_mat[i, j] = d
            dist_mat[j, i] = d
    return dist_mat

def avg_pairwise_distance(dist_mat, idxs):
    """
    Computes average pairwise distance within the subset indexed by idxs
    """
    if len(idxs) < 2:
        return 0.0
    submat = dist_mat[np.ix_(idxs, idxs)]
    tril = submat[np.tril_indices(len(idxs), k=-1)]
    return np.mean(tril)

def test_statistic_indices(dist_mat, idx_group1, idx_group2):
    """
    Robinson & Turner test statistic based on indices and precomputed distance matrix
    """
    return avg_pairwise_distance(dist_mat, idx_group1) + avg_pairwise_distance(dist_mat, idx_group2)

def permutation_test(pd_group1, pd_group2, num_permutations=1000, seed=42):
    """
    Fast permutation test with precomputed Wasserstein distances
    """
    random.seed(seed)
    np.random.seed(seed)

    all_diagrams = pd_group1 + pd_group2
    n1 = len(pd_group1)
    n_total = len(all_diagrams)

    dist_mat = compute_distance_matrix(all_diagrams)

    idxs = np.arange(n_total)
    idx_g1 = np.arange(n1)
    idx_g2 = np.arange(n1, n_total)

    T_obs = test_statistic_indices(dist_mat, idx_g1, idx_g2)

    T_perm = []
    for _ in range(num_permutations):
        np.random.shuffle(idxs)
        new_g1 = idxs[:n1]
        new_g2 = idxs[n1:]
        t = test_statistic_indices(dist_mat, new_g1, new_g2)
        T_perm.append(t)

    T_perm = np.array(T_perm)
    p_value = (np.sum(T_perm <= T_obs) + 1) / (num_permutations + 1)
    return T_obs, T_perm, p_value

""""
================================================================================================================================================================
Persistence Landscape test

The base code was referenced from 
https://persim.scikit-tda.org/en/latest/notebooks/Differentiation%20with%20Persistence%20Landscapes.html#Establish-the-baseline
================================================================================================================================================================
"""

def permutation_pl_test(pl_list1,pl_list2,seed=42,num_perms=1000):
    avg_pl_list1 = average_approx(pl_list1)
    avg_pl_list2 = average_approx(pl_list2)
    [avg_pl_list1_snapped, avg_pl_list2_snapped] = snap_pl([avg_pl_list1, avg_pl_list2])
    true_diff_pl = avg_pl_list1_snapped - avg_pl_list2_snapped
    significance = true_diff_pl.p_norm(p=2)

    comb_pl =pl_list1  + pl_list2
    sig_count = 0
    random.seed(seed)
    for shuffle in range(num_perms):
        A_indices = random.sample(range(len(pl_list1)+len(pl_list2)),len(pl_list1))
        B_indices = [_ for _ in range(len(pl_list1)+len(pl_list2)) if _ not in A_indices]

        A_pl = [comb_pl[i] for i in A_indices]
        B_pl = [comb_pl[j] for j in B_indices]

        A_avg = average_approx(A_pl)
        B_avg = average_approx(B_pl)
        [A_avg_sn, B_avg_sn] = snap_pl([A_avg,B_avg])

        shuff_diff = A_avg_sn - B_avg_sn
        if (shuff_diff.p_norm(p=2) >= significance): sig_count += 1

    pval = (sig_count+1)/(num_perms+1)
    return pval


def landscape_l2_gram(landscapes, start=None, stop=None):
    """Compute the L2 Gram matrix of fixed-grid persistence landscapes.

    Each item in ``landscapes`` may be either a ``PersLandscapeApprox``
    object or a two-dimensional array with shape ``(depth, num_steps)``.
    Landscapes with different depths are padded with zero layers, as in the
    usual persistence-landscape arithmetic.

    All landscapes must already use the same grid.  For
    ``PersLandscapeApprox`` inputs, ``start`` and ``stop`` are inferred from
    the objects and checked for consistency.  For raw arrays, pass ``start``
    and ``stop`` if the absolute L2 inner products are needed.  If they are
    omitted, unit grid spacing is used; this common rescaling does not change
    the permutation p-value returned by ``permutation_pl_test_fast``.

    Parameters
    ----------
    landscapes : sequence
        Fixed-grid persistence landscapes.
    start, stop : float or None
        Common grid endpoints.  Supply both or neither.

    Returns
    -------
    np.ndarray
        Symmetric matrix whose ``(i, j)`` entry is the L2 inner product of
        landscapes ``i`` and ``j``.
    """
    landscapes = list(landscapes)
    if len(landscapes) == 0:
        raise ValueError("At least one persistence landscape is required.")
    if (start is None) != (stop is None):
        raise ValueError("start and stop must be supplied together.")

    values_list = []
    object_grids = []
    num_steps = None
    max_depth = 0

    for landscape_index, landscape in enumerate(landscapes):
        try:
            values = np.asarray(
                getattr(landscape, "values", landscape),
                dtype=float,
            )
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Landscape {landscape_index} does not contain a numeric "
                "fixed-grid value array."
            ) from exc

        if values.ndim != 2:
            raise ValueError(
                f"Landscape {landscape_index} has shape {values.shape}; "
                "expected (depth, num_steps)."
            )
        if values.shape[0] < 1 or values.shape[1] < 2:
            raise ValueError(
                f"Landscape {landscape_index} has shape {values.shape}; "
                "depth must be positive and num_steps must be at least 2."
            )
        if not np.all(np.isfinite(values)):
            raise ValueError(
                f"Landscape {landscape_index} contains non-finite values."
            )

        if num_steps is None:
            num_steps = values.shape[1]
        elif values.shape[1] != num_steps:
            raise ValueError(
                "All landscapes must use the same number of grid points; "
                f"landscape 0 has {num_steps}, whereas landscape "
                f"{landscape_index} has {values.shape[1]}."
            )

        declared_num_steps = getattr(landscape, "num_steps", None)
        if (
            declared_num_steps is not None
            and int(declared_num_steps) != values.shape[1]
        ):
            raise ValueError(
                f"Landscape {landscape_index} declares "
                f"num_steps={declared_num_steps}, but its values have "
                f"{values.shape[1]} grid points."
            )

        landscape_start = getattr(landscape, "start", None)
        landscape_stop = getattr(landscape, "stop", None)
        if (landscape_start is None) != (landscape_stop is None):
            raise ValueError(
                f"Landscape {landscape_index} has incomplete grid metadata."
            )
        if landscape_start is not None:
            object_grids.append(
                (landscape_index, float(landscape_start), float(landscape_stop))
            )

        values_list.append(np.ascontiguousarray(values, dtype=float))
        max_depth = max(max_depth, values.shape[0])

    if start is None:
        if object_grids:
            _, grid_start, grid_stop = object_grids[0]
        else:
            # Raw arrays do not carry their physical grid.  Unit spacing is
            # sufficient for the permutation p-value because every statistic
            # is multiplied by the same positive constant.
            grid_start = 0.0
            grid_stop = float(num_steps - 1)
    else:
        grid_start = float(start)
        grid_stop = float(stop)

    if not np.isfinite(grid_start) or not np.isfinite(grid_stop):
        raise ValueError("start and stop must be finite.")
    if grid_stop <= grid_start:
        raise ValueError("stop must be greater than start.")

    grid_tolerance = (
        64.0
        * np.finfo(float).eps
        * max(1.0, abs(grid_start), abs(grid_stop))
    )
    for landscape_index, object_start, object_stop in object_grids:
        if (
            abs(object_start - grid_start) > grid_tolerance
            or abs(object_stop - grid_stop) > grid_tolerance
        ):
            raise ValueError(
                "All landscapes must use one common grid; landscape "
                f"{landscape_index} uses [{object_start}, {object_stop}], "
                f"not [{grid_start}, {grid_stop}]."
            )

    number_of_landscapes = len(values_list)
    interval_count = num_steps - 1
    feature_count = max_depth * interval_count

    left_endpoints = np.zeros(
        (number_of_landscapes, feature_count),
        dtype=float,
    )
    right_endpoints = np.zeros_like(left_endpoints)

    for landscape_index, values in enumerate(values_list):
        used = values.shape[0] * interval_count
        left_endpoints[landscape_index, :used] = values[:, :-1].reshape(-1)
        right_endpoints[landscape_index, :used] = values[:, 1:].reshape(-1)

    grid_step = (grid_stop - grid_start) / interval_count
    left_left = left_endpoints @ left_endpoints.T
    right_right = right_endpoints @ right_endpoints.T
    left_right = left_endpoints @ right_endpoints.T

    gram = (grid_step / 6.0) * (
        2.0 * left_left
        + 2.0 * right_right
        + left_right
        + left_right.T
    )

    # Remove harmless BLAS-level asymmetry before quadratic forms.
    return 0.5 * (gram + gram.T)


def permutation_pl_test_fast(
    pl_list1,
    pl_list2,
    seed=42,
    num_perms=1000,
    start=None,
    stop=None,
):
    """Fast permutation test for equality of mean persistence landscapes.

    This is the fixed-grid Gram-matrix version of ``permutation_pl_test``.
    It computes all pairwise L2 inner products once, then evaluates the
    observed and permuted squared distances between sample means as quadratic
    forms.  Thus it avoids repeated averaging and grid snapping inside the
    permutation loop.

    Parameters
    ----------
    pl_list1, pl_list2 : sequence
        Nonempty samples of ``PersLandscapeApprox`` objects or numeric arrays
        with shape ``(depth, num_steps)``.  Every item must use the same grid.
    seed : int or None, default=42
        Seed for a local Python random-number generator.  The global random
        state is not modified.
    num_perms : int, default=1000
        Number of random label permutations.
    start, stop : float or None
        Optional common grid endpoints.  These are normally inferred from
        ``PersLandscapeApprox`` inputs.  They may be omitted for raw arrays
        because their common scale does not affect the p-value.

    Returns
    -------
    float
        Monte Carlo permutation p-value with the plus-one correction.
    """
    pl_list1 = list(pl_list1)
    pl_list2 = list(pl_list2)
    n = len(pl_list1)
    m = len(pl_list2)

    if n == 0 or m == 0:
        raise ValueError("Both PL samples must be nonempty.")
    if isinstance(num_perms, (bool, np.bool_)):
        raise ValueError("num_perms must be a positive integer.")
    try:
        integer_num_perms = int(num_perms)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("num_perms must be a positive integer.") from exc
    if integer_num_perms != num_perms or integer_num_perms < 1:
        raise ValueError("num_perms must be a positive integer.")
    num_perms = integer_num_perms

    combined = pl_list1 + pl_list2
    total = n + m
    gram = landscape_l2_gram(combined, start=start, stop=stop)

    # Row 0 is the observed split; rows 1:B are random splits.
    coefficients = np.full(
        (num_perms + 1, total),
        -1.0 / m,
        dtype=float,
    )
    coefficients[0, :n] = 1.0 / n

    rng = random.Random(seed)
    population = range(total)
    for permutation_index in range(1, num_perms + 1):
        group_x = rng.sample(population, n)
        coefficients[permutation_index, group_x] = 1.0 / n

    squared_statistics = np.einsum(
        "bi,ij,bj->b",
        coefficients,
        gram,
        coefficients,
        optimize=True,
    )
    squared_statistics = np.maximum(squared_statistics, 0.0)

    observed_squared = float(squared_statistics[0])
    numerical_tolerance = (
        64.0
        * np.finfo(float).eps
        * max(1.0, float(np.max(squared_statistics)))
    )
    exceedances = np.count_nonzero(
        squared_statistics[1:]
        >= observed_squared - numerical_tolerance
    )

    return float((exceedances + 1) / (num_perms + 1))

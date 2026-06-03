import math
import os
import time

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("Python benchmark requires numpy") from exc


def simulated_community(nsites=400, ntaxa=200, total=10_000):
    matrix = np.empty((nsites, ntaxa), dtype=float)
    for site in range(1, nsites + 1):
        weights = np.array(
            [1 / (1 + ((taxon + 3 * site) % ntaxa)) ** 1.15 for taxon in range(1, ntaxa + 1)],
            dtype=float,
        )
        values = np.rint(total * weights / weights.sum())
        values[site % ntaxa] = max(values[site % ntaxa], 1)
        matrix[site - 1, :] = values
    return matrix


def best_time(func, repeats=5, inner=1):
    func()
    best = math.inf
    for _ in range(repeats):
        start = time.perf_counter()
        for _ in range(inner):
            func()
        best = min(best, (time.perf_counter() - start) / inner)
    return best


def alpha_diversity_numpy(matrix):
    totals = matrix.sum(axis=1)
    probabilities = matrix / totals[:, None]
    positive = probabilities > 0
    richness = positive.sum(axis=1)
    logs = np.zeros_like(probabilities)
    np.log2(probabilities, out=logs, where=positive)
    shannon = -(np.where(positive, probabilities * logs, 0)).sum(axis=1)
    simpson = (probabilities * probabilities).sum(axis=1)
    return richness, shannon, 2**shannon, simpson


def richness_numpy(matrix):
    return (matrix > 0).sum(axis=1)


def shannon_numpy(matrix):
    totals = matrix.sum(axis=1)
    probabilities = matrix / totals[:, None]
    positive = probabilities > 0
    logs = np.zeros_like(probabilities)
    np.log2(probabilities, out=logs, where=positive)
    return -(np.where(positive, probabilities * logs, 0)).sum(axis=1)


def bray_curtis_scipy(matrix):
    from scipy.spatial.distance import pdist, squareform

    return squareform(pdist(matrix, metric="braycurtis"))


def jaccard_scipy(matrix):
    from scipy.spatial.distance import pdist, squareform

    return squareform(pdist(matrix > 0, metric="jaccard"))


def hellinger_numpy(matrix):
    probabilities = matrix / matrix.sum(axis=1)[:, None]
    roots = np.sqrt(probabilities)
    nsites = roots.shape[0]
    result = np.empty((nsites, nsites), dtype=float)
    for i in range(nsites):
        result[i, i] = 0
        for j in range(i + 1, nsites):
            value = math.sqrt(((roots[i] - roots[j]) ** 2).sum() / 2)
            result[i, j] = value
            result[j, i] = value
    return result


def main():
    nsites = int(os.getenv("DIVERSITY_BENCH_NSITES", "400"))
    ntaxa = int(os.getenv("DIVERSITY_BENCH_NTAXA", "200"))
    total = int(os.getenv("DIVERSITY_BENCH_TOTAL", "10000"))
    repeats = int(os.getenv("DIVERSITY_BENCH_REPEATS", "5"))
    inner = int(os.getenv("DIVERSITY_BENCH_INNER", "1"))
    matrix = simulated_community(nsites, ntaxa, total)

    print("language,package,task,nsites,ntaxa,repeats,inner,best_seconds")
    print(f"Python,numpy,richness,{nsites},{ntaxa},{repeats},{inner},{best_time(lambda: richness_numpy(matrix), repeats, inner)}")
    print(f"Python,numpy,shannon_entropy,{nsites},{ntaxa},{repeats},{inner},{best_time(lambda: shannon_numpy(matrix), repeats, inner)}")
    print(f"Python,numpy,alpha_diversity,{nsites},{ntaxa},{repeats},{inner},{best_time(lambda: alpha_diversity_numpy(matrix), repeats, inner)}")
    try:
        scipy_time = best_time(lambda: bray_curtis_scipy(matrix), repeats, inner)
        print(f"Python,scipy,bray_curtis_distance_matrix,{nsites},{ntaxa},{repeats},{inner},{scipy_time}")
        jaccard_time = best_time(lambda: jaccard_scipy(matrix), repeats, inner)
        print(f"Python,scipy,jaccard_distance_matrix,{nsites},{ntaxa},{repeats},{inner},{jaccard_time}")
    except ImportError:
        print(f"Python,scipy,bray_curtis_distance_matrix,{nsites},{ntaxa},{repeats},{inner},unavailable")
        print(f"Python,scipy,jaccard_distance_matrix,{nsites},{ntaxa},{repeats},{inner},unavailable")
    print(f"Python,numpy,hellinger_distance_matrix,{nsites},{ntaxa},{repeats},{inner},{best_time(lambda: hellinger_numpy(matrix), repeats, inner)}")


if __name__ == "__main__":
    main()

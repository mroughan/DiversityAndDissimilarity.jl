# Benchmarks

The `benchmark/` directory contains small cross-language scripts for comparing
common diversity workflows on a deterministic simulated community matrix. The
matrix is generated at runtime, so no large benchmark dataset is committed and
Git LFS is not needed.

The default benchmark uses 400 sites, 200 taxa, and a total abundance near
10,000 per site. This is large enough for pairwise distance matrices to do real
work while still being practical on a laptop.

## Running The Julia Benchmark

From the repository root:

```bash
julia --project=. benchmark/julia_benchmark.jl
```

The script reports CSV rows for:

- `alpha_diversity`
- `bray_curtis_distance_matrix`
- `hellinger_distance_matrix`

## Running Python And R Comparisons

The Python script requires `numpy`; Bray-Curtis comparison additionally uses
`scipy` when available.

```bash
python3 benchmark/python_benchmark.py
```

The R script requires `Rscript` and the `vegan` package:

```bash
Rscript benchmark/r_vegan_benchmark.R
```

## Generating A Multi-Scenario Report

Run all configured scenarios and write a markdown report plus raw CSV:

```bash
julia --project=. benchmark/generate_report.jl
```

The report generator uses `~/.venvs/diversity-bench/bin/python` when present,
or `python3` otherwise. Set `DIVERSITY_BENCH_PYTHON` to override this.
The generated report includes processor, platform, Julia version, Julia thread
count, and Python interpreter details so timing results can be interpreted in
context.

The default scenarios include a small smoke-test matrix, the default
400-site/200-taxon matrix, and larger many-site, wide, and large matrices.
Generated output is written to:

- `benchmark/results/benchmark-results.csv`
- `benchmark/results/benchmark-report.md`

All scripts accept the same environment variables:

```bash
DIVERSITY_BENCH_NSITES=800 \
DIVERSITY_BENCH_NTAXA=300 \
DIVERSITY_BENCH_TOTAL=10000 \
DIVERSITY_BENCH_REPEATS=3 \
DIVERSITY_BENCH_INNER=1 \
julia --project=. benchmark/julia_benchmark.jl
```

`DIVERSITY_BENCH_INNER` controls how many times each timed expression is run
inside one measurement. The scripts report per-call time. This is useful for
small scenarios where operating-system or language timer resolution can
otherwise produce zero or noisy elapsed times.

The Julia benchmark runs every timed function once before measurement so JIT
compilation work is not included in reported timings.

## Safe and Pre-validated Pathways

By default every function validates its input on every call — checking for
negative abundances, non-finite values, and all-zero rows before any
computation begins. This is the safe, user-facing default.

When the same community matrix is passed to many functions in a pipeline the
validation overhead accumulates. `validate` runs all checks once and returns a
[`Validated`](@ref) wrapper that signals the data is already clean:

```julia
v = validate(community)   # validates once; returns Validated{Matrix{Float64}}

richness(v)               # computation only — no re-validation
shannon_entropy(v)
alpha_diversity(v)
bray_curtis_distance(v)
jaccard_distance(v)
hellinger_distance(v)
```

The Julia benchmark outputs both pathways as separate CSV rows. Tasks ending in
`_prevalidated` use pre-validated data and are the fair point of comparison
with Python/NumPy and R/vegan, which do not re-validate inputs on each call.
The plain task names include per-call validation and reflect the full default
user experience.

!!! warning "Risk of pre-validated pathway"
    Constructing `Validated(data)` directly — without calling `validate` —
    bypasses all checks. Invalid data may produce silently wrong results.
    Always call `validate(data)` to obtain a `Validated` wrapper from
    unverified input.

## Notes On Interpretation

These benchmarks are intended as workflow checks, not universal language
shootouts. They compare practical calls users are likely to make:
row-wise alpha diversity and dense pairwise distance matrices. Dense pairwise
matrices scale as ``O(n_{sites}^2 n_{taxa})`` in time and
``O(n_{sites}^2)`` in memory regardless of language.

If a benchmark output file becomes useful for releases, generate it from the
scripts and commit only a small summary. Store large raw timing outputs outside
the repository or in Git LFS.

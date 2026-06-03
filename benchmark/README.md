# Benchmark Setup

The benchmark scripts generate a deterministic simulated community matrix at
runtime. No large input data is committed to the repository.

## Julia

```bash
julia --project=. benchmark/julia_benchmark.jl
```

## Python

Create a dedicated environment outside the repository:

```bash
python3 -m venv ~/.venvs/diversity-bench
~/.venvs/diversity-bench/bin/pip install -r benchmark/requirements.txt
~/.venvs/diversity-bench/bin/python benchmark/python_benchmark.py
```

## R

Install the benchmark packages once:

```bash
Rscript benchmark/r-packages.R
```

Then run:

```bash
Rscript benchmark/r_vegan_benchmark.R
```

## Scale

All scripts use the same optional environment variables:

```bash
DIVERSITY_BENCH_NSITES=800 \
DIVERSITY_BENCH_NTAXA=300 \
DIVERSITY_BENCH_TOTAL=10000 \
DIVERSITY_BENCH_REPEATS=3 \
DIVERSITY_BENCH_INNER=1 \
julia --project=. benchmark/julia_benchmark.jl
```

`DIVERSITY_BENCH_INNER` repeats the timed expression inside each measurement
and reports per-call time. Increase it for very small problems to avoid coarse
timer-resolution artifacts.

The Julia benchmark runs every timed function once before measurement so JIT
compilation work is not included in reported timings.

## Full Report

Run the multi-scenario report generator from the repository root:

```bash
julia --project=. benchmark/generate_report.jl
```

By default it uses `~/.venvs/diversity-bench/bin/python` when that environment
exists. Override the Python interpreter with:

```bash
DIVERSITY_BENCH_PYTHON=/path/to/python julia --project=. benchmark/generate_report.jl
```

The generated files are:

- `benchmark/results/benchmark-results.csv`
- `benchmark/results/benchmark-report.md`

The markdown report includes processor, platform, Julia version, Julia thread
count, and Python interpreter details.

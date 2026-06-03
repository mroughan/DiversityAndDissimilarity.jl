# Release Readiness

This page collects the practical checks to run before tagging or registering a
release.

## Metadata

- `Project.toml` has a stable UUID, version, author, description, dependency
  bounds, and Julia compatibility bounds.
- `LICENSE`, `README.md`, `CHANGELOG.md`, and documentation are present.
- Public API additions are exported, documented, and tested.

## Checks

Run the package tests:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Build the documentation:

```bash
julia --project=docs docs/make.jl
```

Run the benchmark report when performance-sensitive code changes:

```bash
benchmark/run_report.sh
```

## Registration

Before registering in Julia General, make sure CI, docs, and quality checks are
green on the release commit. Use a patch/minor version number that matches the
scope of the public API changes.

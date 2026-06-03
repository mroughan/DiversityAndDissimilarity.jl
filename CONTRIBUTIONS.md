# Contributions

Thank you for improving `DiversityAndDissimilarity.jl`. This project is small
on purpose, so the most valuable contributions are clear, convention-aware,
well-tested changes that keep the package easy to install and easy to reason
about.

Before making a substantial change, read `ARCHITECTURE.md`. It describes the
goals, non-goals, framework structure, and design decisions that should be
preserved as the package grows.

## Development Setup

From the repository root:

```bash
julia --project=.
```

Then in the Julia package prompt:

```julia
pkg> instantiate
pkg> test
```

The package supports Julia 1.10 and newer 1.x releases. Runtime dependencies
should remain minimal; optional tooling belongs in scoped environments such as
`docs` or `test/quality`.

## Contribution Workflow

1. Keep changes focused. Separate unrelated refactors, new indices, docs
   updates, and validation data when possible.
2. Prefer the existing dispatch style over string-based method selection.
3. Add or update tests with the code change.
4. Update documentation in the same change when public behavior changes.
5. Run the relevant verification commands before opening a PR.

## Adding Or Changing An Index

New indices should participate in the whole framework, not only expose a
numeric function.

Checklist:

1. Add the index type to `src/diversity.jl` or `src/similarity.jl`.
2. Implement the generic operation first: `entropy`, `diversity`,
   `effective_diversity`, `similarity`, or `dissimilarity`.
3. Add convenience functions only after the generic method is in place.
4. Add metadata in `src/framework.jl`, including family, input/output mode,
   range, interpreted bounds, metric-like traits, formula, aliases, and notes.
5. Add tests for formulas, edge cases, metadata, matrix behavior, and
   convenience aliases.
6. Add documentation in the relevant manual page and API reference.
7. Add validation references when an external package, paper, or textbook gives
   comparable expected values.
8. Regenerate generated documentation assets if the type tree changes:

```bash
julia --project=docs docs/make_type_trees.jl
```

## Metadata Expectations

Every public index should expose meaningful metadata where the answer is known:

- `index_family`
- `input_mode`
- `output_mode`
- `index_range`
- `index_bounds`
- `is_finite`
- `is_symmetric`
- `is_nonnegative`
- `is_bounded`
- `is_metric`
- `is_triangular`
- metric-like descriptors such as `is_semimetric` and `is_premetric`
- `is_similarity`
- `is_dissimilarity`
- `requires_probabilities`
- `supports_matrix_kernel`

Use `:unknown` for unencoded mathematical properties rather than overstating a
claim. This is part of the package's public contract.

## Tests

Run the main test suite:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Run coverage when changing framework or shared behavior:

```bash
julia --project=. -e 'using Pkg; Pkg.test(; coverage=true)'
```

Run quality checks separately:

```bash
julia --project=test/quality -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate(); include("test/quality/runtests.jl")'
```

The main tests should stay independent of optional R, Python, LaTeX, and
Graphviz tooling. Use generated fixtures in `validation` for cross-package
values rather than requiring external packages in normal CI.

## Documentation

Build the manual with:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate(); include("docs/make.jl")'
```

When public behavior changes, update the usual places:

- source docstrings
- README when the change affects the quick-start or headline feature set
- relevant pages in `docs/src`
- `docs/src/api.md`
- `docs/src/framework.md` for trait or metadata changes
- `docs/src/type-structure.md` and generated assets when public types change
- `notes` when methodology or references change

The type-tree diagram is generated from exported package types:

```bash
julia --project=docs docs/make_type_trees.jl
```

This command writes DOT, SVG, and PDF files under `docs/src/assets`.
Graphviz's `dot` command is required for SVG/PDF rendering.

## Validation Data

Use `validation` for external reference datasets, manifests, and generated
expected values. Prefer small fixtures that are:

- openly available,
- stable across upstream package versions,
- easy to cite,
- useful for checking convention-sensitive behavior.

Do not add large datasets unless they are clearly necessary. When possible,
store the raw fixture, the script or manifest that generated expected values,
and the upstream package/version information.

## Style

- Keep code simple and explicit.
- Prefer clear loops and helper functions over clever abstractions.
- Avoid adding dependencies for convenience alone.
- Use structured input handling rather than ad hoc string parsing.
- Keep comments short and reserved for non-obvious logic.
- Preserve ASCII in source and docs unless there is a clear reason otherwise.
- Keep matrix semantics consistent: samples in rows, taxa/categories in
  columns.

## Pull Request Checklist

Before submitting:

- [ ] Main tests pass.
- [ ] Quality tests pass when the change affects public API, exports, or shared
      behavior.
- [ ] Docs build when docstrings or manual pages change.
- [ ] Generated type-tree assets are updated when exported public types change.
- [ ] Metadata and bounds are updated for any new index.
- [ ] Validation references are added or updated when conventions depend on an
      external source.
- [ ] `git diff --check` is clean.

## Questions And Design Changes

For changes that alter conventions, public API shape, dependency policy, or the
meaning of an existing index, open an issue or design discussion first. The
answer should be reflected in `ARCHITECTURE.md` or the relevant documentation
so future contributors can preserve the decision.

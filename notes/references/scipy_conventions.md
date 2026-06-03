# SciPy Convention Comparisons

Reference: SciPy documentation, `scipy.spatial.distance`
https://docs.scipy.org/doc/scipy/reference/spatial.distance.html

This document records where `DiversityAndDissimilarity.jl` agrees with SciPy
and where conventions differ, so tests can document both.

---

## Bray-Curtis: identical formula

Both packages use:
```
BC(u, v) = Σ |uᵢ − vᵢ| / Σ (uᵢ + vᵢ)
```

SciPy reference values (from documentation examples):
```
braycurtis([1,0,0], [0,1,0]) = 1.0          (disjoint support)
braycurtis([1,1,0], [0,1,0]) = 1/3 ≈ 0.3333 (one shared, one extra)
```

Both values are reproduced exactly by this package.

---

## Canberra: convention difference (averaging)

SciPy uses the **unaveraged** sum:
```
d(u, v) = Σ |uᵢ − vᵢ| / (|uᵢ| + |vᵢ|)
```

This package uses the **averaged** form (divided by the number of nonzero terms m):
```
C(x, y) = (1/m) Σᵢ:xᵢ+yᵢ>0  |xᵢ − yᵢ| / (xᵢ + yᵢ)
```

Reference example — [1, 2, 3] vs [2, 1, 0]:

| Term | |u−v| / (|u|+|v|) |
|------|-------------------|
| (1,2) | 1/3               |
| (2,1) | 1/3               |
| (3,0) | 3/3 = 1           |

Sum = 1/3 + 1/3 + 1 = 5/3

- **SciPy**:      5/3 ≈ 1.6667
- **This package**: (5/3)/3 = 5/9 ≈ 0.5556  (averaged over m=3 nonzero terms)

Ratio: SciPy / ours = 3 (for this example).

In general, `scipy_canberra = ours * m` where m is the number of coordinate pairs
where xᵢ + yᵢ > 0.

---

## Jensen-Shannon distance: different logarithm base

SciPy uses the **natural logarithm** and returns √(JSD_nats):
```python
scipy.spatial.distance.jensenshannon(p, q)  # base=e by default
```

This package uses **log base 2** (bits) and returns √(JSD_bits):
```julia
jensen_shannon_distance(p, q)  # base=2 by default
```

For any pair of distributions: `ours = scipy × √(log₂ e) ≈ scipy × 1.2011`

Concretely: `√(JSD_bits) = √(JSD_nats / ln 2) = √(JSD_nats) × √(1/ln 2)`

### Reference examples

| p               | q               | This package (base=2) | SciPy (nats)    |
|-----------------|-----------------|-----------------------|-----------------|
| [1.0, 0.0]      | [0.5, 0.5]      | 0.55792               | 0.46450         |
| [1.0, 0.0, 0.0] | [0.0, 1.0, 0.0] | 1.00000               | 0.83255         |
| [1.0, 0.0, 0.0] | [1.0, 0.0, 0.0] | 0.00000               | 0.00000         |

The ratio for the first two entries: 0.55792/0.46450 ≈ 0.83255⁻¹ ≈ √(1/ln 2) ≈ 1.2011 ✓

SciPy values (0.46450 and 0.83255) are documented in the official SciPy API reference for
`scipy.spatial.distance.jensenshannon`.

Passing `base=ℯ` to this package gives the natural-log form:
```julia
jensen_shannon_distance(p, q; base=ℯ)  # matches scipy
```

---

## Hellinger: consistent (same formula, different wrapper path)

SciPy has no direct `hellinger` in `scipy.spatial.distance`; the Hellinger distance
is obtained via `numpy` as `sqrt(sum((sqrt(p) - sqrt(q))^2) / 2)`.

This package computes the same formula directly. The Chord distance is:
```
chord(p, q) = sqrt(sum((sqrt(p) - sqrt(q))^2))  =  sqrt(2) × hellinger(p, q)
```

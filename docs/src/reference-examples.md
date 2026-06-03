# Reference Examples

This page collects small examples whose expected values match conventions used
by commonly used ecology and biodiversity packages. They are also covered by the
test suite so that shared indices remain comparable over time.

## scikit-bio Community Matrix

The scikit-bio diversity tutorial uses the following six-sample community
matrix. Each row is a sample and each column is a taxon count.

```jldoctest
julia> using DiversityAndDissimilarity

julia> data = [
           [23, 64, 14, 0, 0, 3, 1],
           [0, 3, 35, 42, 0, 12, 1],
           [0, 5, 5, 0, 40, 40, 0],
           [44, 35, 9, 0, 1, 0, 0],
           [0, 2, 8, 0, 35, 45, 1],
           [0, 0, 25, 35, 0, 19, 0],
       ];

julia> richness.(data)
6-element Vector{Int64}:
 5
 5
 4
 4
 5
 3

julia> bray_curtis_distance(data[1], data[2])
0.7878787878787878

julia> bray_curtis_distance(data[1], data[4])
0.30927835051546393

julia> bray_curtis_distance(data[3], data[5])
0.09392265193370165

julia> jaccard_distance(data[1], data[2])
0.33333333333333337

julia> jaccard_similarity(data[1], data[2])
0.6666666666666666
```

## vegan Conventions

The vegan `diversity` documentation defines Shannon entropy with a configurable
logarithm base,

```math
H = -\sum_i p_i \log_b(p_i),
```

and its Simpson-family outputs as

```math
D = \sum_i p_i^2,\qquad
1-D,\qquad
\frac{1}{D}.
```

For `x = [1, 1, 2]`, the relative abundances are
``p = [1/4, 1/4, 1/2]``:

```jldoctest
julia> using DiversityAndDissimilarity

julia> x = [1, 1, 2];

julia> entropy(Shannon(; base=ℯ), x)
1.0397207708399179

julia> gini_simpson_index(x)
0.625

julia> inverse_simpson_index(x)
2.6666666666666665

julia> pielou_evenness(x)
0.946394630357186

julia> fisher_alpha(x) * log1p(sum(x) / fisher_alpha(x))
3.0
```

The vegan `vegdist` documentation describes Bray-Curtis dissimilarity as

```math
\frac{\sum_i |x_i - y_i|}{\sum_i (x_i + y_i)}.
```

```jldoctest
julia> using DiversityAndDissimilarity

julia> bray_curtis_distance([1, 2, 3], [2, 2, 0])
0.4
```

Presence/absence pairwise conventions:

```jldoctest
julia> using DiversityAndDissimilarity

julia> left = [1, 1, 0, 1];

julia> right = [1, 0, 1, 1];

julia> jaccard_similarity(left, right)
0.5

julia> sorensen_dice_index(left, right)
0.6666666666666666
```

## Magurran (2004) Five-Species Community

A five-species community used as a worked example in Magurran (2004)
*Measuring Biological Diversity* (Blackwell). The relative abundances are
``p = [135, 76, 45, 22, 13] / 291``.

```jldoctest
julia> using DiversityAndDissimilarity

julia> x = [135, 76, 45, 22, 13];

julia> richness(x)
5

julia> simpson_index(x) ≈ (135^2 + 76^2 + 45^2 + 22^2 + 13^2) / 291^2
true

julia> entropy(Shannon(; base=ℯ), x)
1.3296981241766592

julia> pielou_evenness(x) ≈ entropy(Shannon(; base=ℯ), x) / log(5)
true
```

Verify in R with `vegan::diversity(c(135,76,45,22,13), "shannon")`.

## iNEXT Spider Dataset

Spider pitfall-trap data from two hemlock-forest canopy treatments at Harvard
Forest, originally reported by Sackett et al. (2011) and distributed with the
iNEXT package by Chao et al. (2014). The vectors below are the exact abundance
counts extracted from the iNEXT package source.

```jldoctest
julia> using DiversityAndDissimilarity

julia> girdled = [46, 22, 17, 15, 15, 9, 8, 6, 6, 4, 2, 2, 2, 2,
                  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

julia> richness(girdled)
26

julia> chao1(girdled)
39.2

julia> sample_coverage(girdled) ≈ 13/14
true

julia> entropy(Shannon(; base=ℯ), girdled)
2.4898652106915704

julia> logged = [88, 22, 16, 15, 13, 10, 8, 8, 7, 7, 7, 5, 4, 4, 4,
                 3, 3, 3, 3, 2, 2, 2, 2,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

julia> richness(logged)
37

julia> chao1(logged)
55.2

julia> sample_coverage(logged) ≈ 17/18
true

julia> entropy(Shannon(; base=ℯ), logged)
2.6686855819621367
```

Verify Chao1 and ACE in R with `estimateR(spider$Girdled)` from the vegan
package. Full documentation is in `notes/references/spider_dataset.md`.

## SciPy Convention Comparisons

This package and `scipy.spatial.distance` agree on Bray-Curtis but differ on
Canberra (averaging) and Jensen-Shannon (logarithm base).

### Bray-Curtis: identical

```jldoctest
julia> using DiversityAndDissimilarity

julia> bray_curtis_dissimilarity([1, 0, 0], [0, 1, 0])
1.0

julia> bray_curtis_dissimilarity([1, 1, 0], [0, 1, 0])
0.3333333333333333
```

### Canberra: this package averages by number of nonzero terms; SciPy does not

```jldoctest
julia> using DiversityAndDissimilarity

julia> canberra_distance([1, 2, 3], [2, 1, 0])
0.5555555555555555
```

The SciPy unaveraged sum is `5/3 ≈ 1.6667`; this package divides by `m = 3`
nonzero terms to give `5/9 ≈ 0.5556`. See Legendre & Legendre (2012)
*Numerical Ecology* §7.4.

### Jensen-Shannon: base-2 (this package) vs natural log (SciPy)

```jldoctest
julia> using DiversityAndDissimilarity

julia> jensen_shannon_distance([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
1.0

julia> jensen_shannon_distance([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]; base=ℯ)
0.8325546111576977
```

The first call uses `base=2` (default): maximum JSD = 1 bit, so the distance
is 1. The second matches SciPy's natural-log convention. To always match SciPy,
pass `base=ℯ`. The ratio between the two conventions is ``\sqrt{\log_2 e}``.

Full documentation is in `notes/references/scipy_conventions.md`.

## Sources

- [scikit-bio diversity tutorial](https://scikit.bio/docs/latest/diversity.html)
- [vegan diversity documentation](https://vegandevs.github.io/vegan/reference/diversity.html)
- [vegan vegdist documentation](https://vegandevs.github.io/vegan/reference/vegdist.html)
- Magurran, A.E. (2004) *Measuring Biological Diversity*. Blackwell Publishing.
- Sackett et al. (2011) Can. J. Forest Res. 41(2):394–409. doi:10.1139/X10-207
- Chao et al. (2014) Ecol. Monogr. 84(1):45–67. doi:10.1890/13-0133.1
- SciPy 1.0 (2020) Nature Methods 17:261–272. doi:10.1038/s41592-019-0686-2

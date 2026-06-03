# Reference Examples

This page collects small examples whose expected values match conventions used
by commonly used ecology and biodiversity packages. They are also covered by the
test suite so that shared indices remain comparable over time.

## scikit-bio Community Matrix

The scikit-bio diversity tutorial uses the following six-sample community
matrix. Each row is a sample and each column is a taxon count.

```jldoctest
julia> using DiversityIndices

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
julia> using DiversityIndices

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
julia> using DiversityIndices

julia> bray_curtis_distance([1, 2, 3], [2, 2, 0])
0.4
```

Presence/absence pairwise conventions:

```jldoctest
julia> using DiversityIndices

julia> left = [1, 1, 0, 1];

julia> right = [1, 0, 1, 1];

julia> jaccard_similarity(left, right)
0.5

julia> sorensen_dice_index(left, right)
0.6666666666666666
```

## Sources

- [scikit-bio diversity tutorial](https://scikit.bio/docs/latest/diversity.html)
- [vegan diversity documentation](https://vegandevs.github.io/vegan/reference/diversity.html)
- [vegan vegdist documentation](https://vegandevs.github.io/vegan/reference/vegdist.html)

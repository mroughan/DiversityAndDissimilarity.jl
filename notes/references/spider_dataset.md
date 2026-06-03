# iNEXT Spider Dataset

## Source

**Original data:** Sackett, T. E., Record, S., Bewick, S., Baiser, B., Sanders, N. J.,
and Ellison, A. M. (2011). "Disturbance determines sessile and mobile invertebrate
and plant species diversity and composition in New England hemlock forests."
*Canadian Journal of Forest Research*, 41(2):394–409.
doi:10.1139/X10-207

**Package distribution:** Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L.,
Ma, K. H., Colwell, R. K., and Ellison, A. M. (2014). "Rarefaction and extrapolation
with Hill numbers: a framework for sampling and estimation in species diversity
studies." *Ecological Monographs*, 84(1):45–67. doi:10.1890/13-0133.1

Included in the iNEXT R package (Johnson Hsieh et al.) as the bundled `spider`
example dataset; exact abundance vectors extracted from the package source at:
https://github.com/JohnsonHsieh/iNEXT

## Experimental context

Two canopy manipulation treatments of hemlock forests at Harvard Forest (Petersham,
Massachusetts, USA): **Girdled** (hemlock trees girdled, causing gradual die-back)
and **Logged** (hemlock clear-cut). Spider abundance was recorded by standardised
pitfall trapping.

## Abundance vectors

### Girdled (26 species, 168 individuals)

```julia
girdled = [46, 22, 17, 15, 15, 9, 8, 6, 6, 4, 2, 2, 2, 2,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
```

### Logged (37 species, 252 individuals)

```julia
logged = [88, 22, 16, 15, 13, 10, 8, 8, 7, 7, 7, 5, 4, 4, 4,
          3, 3, 3, 3, 2, 2, 2, 2,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
```

## Expected values (computed by DiversityAndDissimilarity.jl, verifiable in R)

To reproduce in R:
```r
library(iNEXT); data(spider)
library(vegan)
vegan::diversity(spider$Girdled, index = "shannon")   # natural log
vegan::diversity(spider$Logged,  index = "shannon")
estimateR(spider$Girdled)   # Chao1, ACE
estimateR(spider$Logged)
```

| Statistic                     | Girdled                   | Logged                    |
|-------------------------------|---------------------------|---------------------------|
| Total individuals (n)         | 168                       | 252                       |
| Observed richness (S)         | 26                        | 37                        |
| Singletons (f₁)               | 12                        | 14                        |
| Doubletons (f₂)               | 4                         | 4                         |
| Sample coverage               | 13/14 ≈ 0.928571          | 17/18 ≈ 0.944444          |
| Chao1 (bias-corrected)        | 39.2 (exact)              | 55.2 (exact)              |
| ACE                           | 48.41522903033908         | 52.68499427262313         |
| Magurran x Shannon H (nats)   | n/a                       | n/a                       |
| Shannon H (nats, base=e)      | 2.4898652106915704        | 2.6686855819621367        |
| Shannon H (bits, base=2)      | 3.592116191946683         | 3.850099454788652         |
| Hill number q=1 (base=e)      | 12.05965049893833         | 14.42100150529583         |
| Hill number q=2               | 7.840000000000005         | 6.761499148211247         |
| Pielou evenness J             | 0.7642085437417918        | 0.7390601632391848        |

### Derivation of exact Chao1 values

Chao1 (bias-corrected) = S_obs + f₁(f₁−1) / (2(f₂+1))

Girdled: 26 + 12×11 / (2×5) = 26 + 132/10 = **39.2** (exact)

Logged:  37 + 14×13 / (2×5) = 37 + 182/10 = **55.2** (exact)

### Derivation of exact sample coverage

Sample coverage = 1 − f₁/n (Good–Turing estimator, Good 1953)

Girdled: 1 − 12/168 = 156/168 = **13/14** (exact)

Logged:  1 − 14/252 = 238/252 = **17/18** (exact)

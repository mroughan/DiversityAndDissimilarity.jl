# DiversityAndDissimilarity Benchmark Report

Generated: 2026-05-23 17:20:16

Python interpreter: `/home/matt/.venvs/diversity-bench/bin/python`

## Environment

| Item | Value |
|---|---|
| CPU | `arrowlake-s` |
| CPU model | `Intel(R) Core(TM) Ultra 9 285K` |
| Machine | `x86_64-linux-gnu` |
| Operating system | `Linux` |
| Julia | `1.12.6` |
| Julia threads | `1` |
| Python interpreter | `/home/matt/.venvs/diversity-bench/bin/python` |

The benchmark matrix is generated deterministically at runtime for each scenario. Timings are best per-call elapsed time over the listed repeats and inner repetitions.

## Scenarios

| Scenario | Sites | Taxa | Total abundance | Repeats | Inner repetitions |
|---|---:|---:|---:|---:|---:|
| small | 120 | 80 | 5000 | 5 | 100 |
| default | 400 | 200 | 10000 | 5 | 20 |
| many_sites | 800 | 250 | 10000 | 3 | 3 |
| wide | 300 | 1000 | 20000 | 3 | 3 |
| large | 1200 | 300 | 20000 | 2 | 1 |

## Figures

The figures use a log-scaled time axis so small row-wise summaries and larger pairwise matrix calculations can be read on the same page. The raw values remain in the results table below.

### Richness Across Scenarios

![Richness Across Scenarios](figures/richness.svg)

Row-wise observed richness.

### Shannon Entropy Across Scenarios

![Shannon Entropy Across Scenarios](figures/shannon_entropy.svg)

Row-wise Shannon entropy.

### Alpha-Diversity Summary Across Scenarios

![Alpha-Diversity Summary Across Scenarios](figures/alpha_diversity.svg)

Compact exploratory alpha-diversity summary. R/vegan is not included because the R benchmark uses separate vegan calls.

### Bray-Curtis Distance Matrix Across Scenarios

![Bray-Curtis Distance Matrix Across Scenarios](figures/bray_curtis_distance_matrix.svg)

Dense pairwise Bray-Curtis dissimilarity matrix.

### Jaccard Distance Matrix Across Scenarios

![Jaccard Distance Matrix Across Scenarios](figures/jaccard_distance_matrix.svg)

Dense pairwise incidence Jaccard distance matrix.

### Hellinger Distance Matrix Across Scenarios

![Hellinger Distance Matrix Across Scenarios](figures/hellinger_distance_matrix.svg)

Dense pairwise Hellinger distance matrix. R/vegan is not included for this direct helper comparison.

## Results

| Scenario | Language | Package | Task | Sites | Taxa | Repeats | Inner | Best seconds |
|---|---|---|---|---:|---:|---:|---:|---:|
| small | Julia | DiversityAndDissimilarity.jl | richness | 120 | 80 | 5 | 100 | 6.33443e-6 |
| small | Julia | DiversityAndDissimilarity.jl | shannon_entropy | 120 | 80 | 5 | 100 | 3.531551e-5 |
| small | Julia | DiversityAndDissimilarity.jl | alpha_diversity | 120 | 80 | 5 | 100 | 0.00024690412 |
| small | Julia | DiversityAndDissimilarity.jl | bray_curtis_distance_matrix | 120 | 80 | 5 | 100 | 0.00015949824 |
| small | Julia | DiversityAndDissimilarity.jl | jaccard_distance_matrix | 120 | 80 | 5 | 100 | 2.399819e-5 |
| small | Julia | DiversityAndDissimilarity.jl | hellinger_distance_matrix | 120 | 80 | 5 | 100 | 0.00018117313 |
| small | Python | numpy | richness | 120 | 80 | 5 | 100 | 5.797159392386675e-06 |
| small | Python | numpy | shannon_entropy | 120 | 80 | 5 | 100 | 7.761111948639154e-05 |
| small | Python | numpy | alpha_diversity | 120 | 80 | 5 | 100 | 8.863512892276048e-05 |
| small | Python | scipy | bray_curtis_distance_matrix | 120 | 80 | 5 | 100 | 0.00022819662000983954 |
| small | Python | scipy | jaccard_distance_matrix | 120 | 80 | 5 | 100 | 0.000533065889030695 |
| small | Python | numpy | hellinger_distance_matrix | 120 | 80 | 5 | 100 | 0.006789145818911493 |
| small | R | vegan | richness | 120 | 80 | 5 | 100 | 7.9999999999999e-05 |
| small | R | vegan | shannon_entropy | 120 | 80 | 5 | 100 | 0.000379999999999998 |
| small | R | vegan | bray_curtis_distance_matrix | 120 | 80 | 5 | 100 | 0.000339999999999998 |
| small | R | vegan | jaccard_distance_matrix | 120 | 80 | 5 | 100 | 0.000509999999999997 |
| default | Julia | DiversityAndDissimilarity.jl | richness | 400 | 200 | 5 | 20 | 5.7699249999999995e-5 |
| default | Julia | DiversityAndDissimilarity.jl | shannon_entropy | 400 | 200 | 5 | 20 | 0.0003003264 |
| default | Julia | DiversityAndDissimilarity.jl | alpha_diversity | 400 | 200 | 5 | 20 | 0.00109443295 |
| default | Julia | DiversityAndDissimilarity.jl | bray_curtis_distance_matrix | 400 | 200 | 5 | 20 | 0.0048645998 |
| default | Julia | DiversityAndDissimilarity.jl | jaccard_distance_matrix | 400 | 200 | 5 | 20 | 0.00031910660000000003 |
| default | Julia | DiversityAndDissimilarity.jl | hellinger_distance_matrix | 400 | 200 | 5 | 20 | 0.0058180327499999995 |
| default | Python | numpy | richness | 400 | 200 | 5 | 20 | 2.9245996847748758e-05 |
| default | Python | numpy | shannon_entropy | 400 | 200 | 5 | 20 | 0.0007863275008276105 |
| default | Python | numpy | alpha_diversity | 400 | 200 | 5 | 20 | 0.0009451454039663077 |
| default | Python | scipy | bray_curtis_distance_matrix | 400 | 200 | 5 | 20 | 0.0060604438534937796 |
| default | Python | scipy | jaccard_distance_matrix | 400 | 200 | 5 | 20 | 0.014508827542886137 |
| default | Python | numpy | hellinger_distance_matrix | 400 | 200 | 5 | 20 | 0.08353081060340628 |
| default | R | vegan | richness | 400 | 200 | 5 | 20 | 4e-04 |
| default | R | vegan | shannon_entropy | 400 | 200 | 5 | 20 | 0.0018 |
| default | R | vegan | bray_curtis_distance_matrix | 400 | 200 | 5 | 20 | 0.00814999999999999 |
| default | R | vegan | jaccard_distance_matrix | 400 | 200 | 5 | 20 | 0.00945 |
| many_sites | Julia | DiversityAndDissimilarity.jl | richness | 800 | 250 | 3 | 3 | 0.00016149233333333333 |
| many_sites | Julia | DiversityAndDissimilarity.jl | shannon_entropy | 800 | 250 | 3 | 3 | 0.0007674520000000001 |
| many_sites | Julia | DiversityAndDissimilarity.jl | alpha_diversity | 800 | 250 | 3 | 3 | 0.0024832976666666665 |
| many_sites | Julia | DiversityAndDissimilarity.jl | bray_curtis_distance_matrix | 800 | 250 | 3 | 3 | 0.025966307666666664 |
| many_sites | Julia | DiversityAndDissimilarity.jl | jaccard_distance_matrix | 800 | 250 | 3 | 3 | 0.0012407106666666667 |
| many_sites | Julia | DiversityAndDissimilarity.jl | hellinger_distance_matrix | 800 | 250 | 3 | 3 | 0.026186475333333334 |
| many_sites | Python | numpy | richness | 800 | 250 | 3 | 3 | 6.377371028065681e-05 |
| many_sites | Python | numpy | shannon_entropy | 800 | 250 | 3 | 3 | 0.0021150296864410243 |
| many_sites | Python | numpy | alpha_diversity | 800 | 250 | 3 | 3 | 0.002373564289882779 |
| many_sites | Python | scipy | bray_curtis_distance_matrix | 800 | 250 | 3 | 3 | 0.030925075989216566 |
| many_sites | Python | scipy | jaccard_distance_matrix | 800 | 250 | 3 | 3 | 0.07298057603960235 |
| many_sites | Python | numpy | hellinger_distance_matrix | 800 | 250 | 3 | 3 | 0.34403428132645786 |
| many_sites | R | vegan | richness | 800 | 250 | 3 | 3 | 0.001 |
| many_sites | R | vegan | shannon_entropy | 800 | 250 | 3 | 3 | 0.005 |
| many_sites | R | vegan | bray_curtis_distance_matrix | 800 | 250 | 3 | 3 | 0.0523333333333333 |
| many_sites | R | vegan | jaccard_distance_matrix | 800 | 250 | 3 | 3 | 0.0563333333333333 |
| wide | Julia | DiversityAndDissimilarity.jl | richness | 300 | 1000 | 3 | 3 | 0.0002336193333333333 |
| wide | Julia | DiversityAndDissimilarity.jl | shannon_entropy | 300 | 1000 | 3 | 3 | 0.0011446146666666667 |
| wide | Julia | DiversityAndDissimilarity.jl | alpha_diversity | 300 | 1000 | 3 | 3 | 0.002219935666666667 |
| wide | Julia | DiversityAndDissimilarity.jl | bray_curtis_distance_matrix | 300 | 1000 | 3 | 3 | 0.01782702833333333 |
| wide | Julia | DiversityAndDissimilarity.jl | jaccard_distance_matrix | 300 | 1000 | 3 | 3 | 0.0007254226666666666 |
| wide | Julia | DiversityAndDissimilarity.jl | hellinger_distance_matrix | 300 | 1000 | 3 | 3 | 0.01803770966666667 |
| wide | Python | numpy | richness | 300 | 1000 | 3 | 3 | 8.398472952346007e-05 |
| wide | Python | numpy | shannon_entropy | 300 | 1000 | 3 | 3 | 0.003316280354435245 |
| wide | Python | numpy | alpha_diversity | 300 | 1000 | 3 | 3 | 0.004256224337344368 |
| wide | Python | scipy | bray_curtis_distance_matrix | 300 | 1000 | 3 | 3 | 0.016909837955608964 |
| wide | Python | scipy | jaccard_distance_matrix | 300 | 1000 | 3 | 3 | 0.04073096734161178 |
| wide | Python | numpy | hellinger_distance_matrix | 300 | 1000 | 3 | 3 | 0.06741864172120889 |
| wide | R | vegan | richness | 300 | 1000 | 3 | 3 | 0.00133333333333333 |
| wide | R | vegan | shannon_entropy | 300 | 1000 | 3 | 3 | 0.00666666666666667 |
| wide | R | vegan | bray_curtis_distance_matrix | 300 | 1000 | 3 | 3 | 0.026 |
| wide | R | vegan | jaccard_distance_matrix | 300 | 1000 | 3 | 3 | 0.0316666666666667 |
| large | Julia | DiversityAndDissimilarity.jl | richness | 1200 | 300 | 2 | 1 | 0.000259373 |
| large | Julia | DiversityAndDissimilarity.jl | shannon_entropy | 1200 | 300 | 2 | 1 | 0.001347677 |
| large | Julia | DiversityAndDissimilarity.jl | alpha_diversity | 1200 | 300 | 2 | 1 | 0.003981658 |
| large | Julia | DiversityAndDissimilarity.jl | bray_curtis_distance_matrix | 1200 | 300 | 2 | 1 | 0.086633236 |
| large | Julia | DiversityAndDissimilarity.jl | jaccard_distance_matrix | 1200 | 300 | 2 | 1 | 0.003480309 |
| large | Julia | DiversityAndDissimilarity.jl | hellinger_distance_matrix | 1200 | 300 | 2 | 1 | 0.072561824 |
| large | Python | numpy | richness | 1200 | 300 | 2 | 1 | 0.0001261490397155285 |
| large | Python | numpy | shannon_entropy | 1200 | 300 | 2 | 1 | 0.004210239974781871 |
| large | Python | numpy | alpha_diversity | 1200 | 300 | 2 | 1 | 0.0046751489862799644 |
| large | Python | scipy | bray_curtis_distance_matrix | 1200 | 300 | 2 | 1 | 0.08314802800305188 |
| large | Python | scipy | jaccard_distance_matrix | 1200 | 300 | 2 | 1 | 0.19697523792274296 |
| large | Python | numpy | hellinger_distance_matrix | 1200 | 300 | 2 | 1 | 0.7872653470840305 |
| large | R | vegan | richness | 1200 | 300 | 2 | 1 | 0.002 |
| large | R | vegan | shannon_entropy | 1200 | 300 | 2 | 1 | 0.00800000000000001 |
| large | R | vegan | bray_curtis_distance_matrix | 1200 | 300 | 2 | 1 | 0.129 |
| large | R | vegan | jaccard_distance_matrix | 1200 | 300 | 2 | 1 | 0.138 |

## Notes

- Dense pairwise distance matrices scale quadratically in the number of sites.
- Reported times are per call. Small scenarios use larger inner repetition counts to avoid coarse timer-resolution artifacts.
- The Julia benchmark runs every timed function once before measurement so JIT compilation work is not included in the reported timings.
- Julia and R timings call package-level APIs. Python uses NumPy/SciPy reference workflows in `benchmark/python_benchmark.py`.
- The benchmark is intended to compare practical workflows, not isolated kernel implementations.

if (!requireNamespace("vegan", quietly = TRUE)) {
  stop("R benchmark requires the vegan package", call. = FALSE)
}

simulated_community <- function(nsites = 400, ntaxa = 200, total = 10000) {
  matrix <- matrix(0, nrow = nsites, ncol = ntaxa)
  for (site in seq_len(nsites)) {
    weights <- vapply(seq_len(ntaxa), function(taxon) {
      1 / (1 + ((taxon + 3 * site) %% ntaxa)) ^ 1.15
    }, numeric(1))
    values <- round(total * weights / sum(weights))
    values[(site %% ntaxa) + 1] <- max(values[(site %% ntaxa) + 1], 1)
    matrix[site, ] <- values
  }
  matrix
}

best_time <- function(expr, repeats = 5, inner = 1) {
  call <- substitute(expr)
  eval(call, parent.frame())
  best <- Inf
  for (i in seq_len(repeats)) {
    elapsed <- system.time({
      for (j in seq_len(inner)) {
        eval(call, parent.frame())
      }
    })[["elapsed"]] / inner
    best <- min(best, elapsed)
  }
  best
}

nsites <- as.integer(Sys.getenv("DIVERSITY_BENCH_NSITES", "400"))
ntaxa <- as.integer(Sys.getenv("DIVERSITY_BENCH_NTAXA", "200"))
total <- as.integer(Sys.getenv("DIVERSITY_BENCH_TOTAL", "10000"))
repeats <- as.integer(Sys.getenv("DIVERSITY_BENCH_REPEATS", "5"))
inner <- as.integer(Sys.getenv("DIVERSITY_BENCH_INNER", "1"))
community <- simulated_community(nsites, ntaxa, total)

cat("language,package,task,nsites,ntaxa,repeats,inner,best_seconds\n")
cat(paste("R", "vegan", "richness", nsites, ntaxa, repeats,
  inner, best_time(vegan::specnumber(community), repeats, inner),
  sep = ","
), "\n")
cat(paste("R", "vegan", "shannon_entropy", nsites, ntaxa, repeats,
  inner, best_time(vegan::diversity(community, index = "shannon", base = 2), repeats, inner),
  sep = ","
), "\n")
cat(paste("R", "vegan", "bray_curtis_distance_matrix", nsites, ntaxa, repeats,
  inner, best_time(vegan::vegdist(community, method = "bray"), repeats, inner),
  sep = ","
), "\n")
cat(paste("R", "vegan", "jaccard_distance_matrix", nsites, ntaxa, repeats,
  inner, best_time(vegan::vegdist(community, method = "jaccard", binary = TRUE), repeats, inner),
  sep = ","
), "\n")

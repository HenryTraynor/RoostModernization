
#' Generate balanced landscapes across Clark-Evans R bins
#'
#' Generates spatial landscapes until each requested Clark-Evans R bin
#' contains n_per_bin landscapes. Each accepted landscape is reused across
#' all values of prop_modern, preserving its coordinates and modernization
#' order across those scenarios.
#'
#' @param n_per_bin Number of base landscapes to generate per R bin.
#' @param R_breaks Numeric vector defining Clark-Evans R bin boundaries.
#' @param prop_modern_values Proportions of modern roosts to evaluate.
#' @param num_roosts Number of roosts per landscape.
#' @param cluster_values Candidate numbers of spatial clusters.
#' @param sd_values Candidate spatial dispersion parameters.
#' @param modern_N_max Capacity assigned to modern roosts.
#' @param traditional_N_max Capacity assigned to traditional roosts.
#' @param min_N_max Minimum permitted roost capacity.
#' @param max_attempts Maximum number of candidate landscapes to generate.
#'
#' @return A list containing landscapes, design, counts, and R_breaks.
#' @export
generate_landscapes <- function(
    n_per_bin,
    R_breaks,
    prop_modern_values = seq(0, 1, by = 0.2),
    num_roosts = 100,
    cluster_values,
    sd_values,
    modern_N_max = 100,
    traditional_N_max = 20,
    min_N_max = 5,
    max_attempts = 10000
) {
  
  n_bins <- length(R_breaks) - 1L
  counts <- integer(n_bins)
  landscapes <- list()
  landscape_index <- 0L
  attempts <- 0L
  
  # Generate candidates until each bin has enough landscapes
  while (any(counts < n_per_bin) && attempts < max_attempts) {
    attempts <- attempts + 1L
    
    num_clusters <- sample(cluster_values, 1)
    sd_value <- sample(sd_values, 1)
    
    # Generate spatial coordinates independently of modernization.
    base_roosts <- create_roosts(
      num_roosts = num_roosts,
      num_clusters = num_clusters,
      sd = sd_value,
      prop_modern = 0,
      min_N_max = min_N_max
    )
    
    R_value <- calculate_clark_evans_R(
      base_roosts$x,
      base_roosts$y
    )
    
    R_bin <- findInterval(
      R_value,
      R_breaks,
      rightmost.closed = TRUE
    )
    
    # Ignore candidates outside the requested R range or in full bins
    if (R_bin < 1 || R_bin > n_bins || counts[R_bin] >= n_per_bin) {
      next
    }
    
    landscape_index <- landscape_index + 1L
    counts[R_bin] <- counts[R_bin] + 1L
    
    # Fix modernization order for this base landscape.
    # The same roosts become modern first at every proportion.
    modernization_order <- sample.int(nrow(base_roosts))
    
    # Generate fixed, normally distributed potential capacities
    # for each roost under each building type.
    
    N_max_traditional <- pmax(
      min_N_max,
      round(rnorm(
        nrow(base_roosts),
        mean = traditional_N_max,
        sd = 0.25 * traditional_N_max
      ))
    )
    
    N_max_modern <- pmax(
      min_N_max,
      round(rnorm(
        nrow(base_roosts),
        mean = modern_N_max,
        sd = 0.25 * modern_N_max
      ))
    )
    
    base_roosts$N_max_traditional <- N_max_traditional
    base_roosts$N_max_modern <- N_max_modern
    base_roosts$modernization_order <- match(
      seq_len(nrow(base_roosts)),
      modernization_order
    )
    
    landscapes[[landscape_index]] <- list(
      landscape_id = landscape_index,
      base_roosts = base_roosts,
      modernization_order = modernization_order,
      Clark_Evans_R = R_value,
      num_clusters = num_clusters,
      sd = sd_value,
      R_bin = R_bin
    )
  }
  
  if (any(counts < n_per_bin)) {
    warning(
      "Could not fill every R bin within max_attempts. ",
      "Try increasing max_attempts or widening R_breaks."
    )
  }
  
  # Cross accepted landscapes with modernization proportions
  design_grid <- expand.grid(
    landscape_index = seq_along(landscapes),
    prop_modern = prop_modern_values,
    KEEP.OUT.ATTRS = FALSE
  )
  
  design <- data.frame(
    landscape_index = design_grid$landscape_index,
    prop_modern = design_grid$prop_modern
  )
  
  design$landscape_id <- NA_integer_
  design$Clark_Evans_R <- NA_real_
  design$num_clusters <- NA_integer_
  design$sd <- NA_real_
  design$R_bin <- NA_integer_
  design$replicate <- NA_integer_
  design$roost_info <- vector("list", nrow(design))
  
  for (i in seq_len(nrow(design))) {
    idx <- design$landscape_index[i]
    landscape <- landscapes[[idx]]
    roost_info <- landscape$base_roosts
    
    n_modern <- round(
      design$prop_modern[i] * nrow(roost_info)
    )
    
    modern_rows <- head(
      landscape$modernization_order,
      n_modern
    )
    
    roost_info$roost_type <- "traditional"
    roost_info$roost_type[modern_rows] <- "modern"
    
    roost_info$N_max <- roost_info$N_max_traditional
    roost_info$N_max[modern_rows] <-
      roost_info$N_max_modern[modern_rows]
    
    # Initialize an uninfected landscape.
    roost_info$S <- 0L
    roost_info$I <- 0L
    roost_info$R <- 0L
    roost_info$N <- 0L
    
    # Remove helper columns from the public roost data.
    roost_info$N_max_traditional <- NULL
    roost_info$N_max_modern <- NULL
    roost_info$modernization_order <- NULL
    
    design$landscape_id[i] <- landscape$landscape_id
    design$Clark_Evans_R[i] <- landscape$Clark_Evans_R
    design$num_clusters[i] <- landscape$num_clusters
    design$sd[i] <- landscape$sd
    design$R_bin[i] <- landscape$R_bin
    design$replicate[i] <- landscape$landscape_id
    design$roost_info[[i]] <- roost_info
  }
  
  list(
    landscapes = landscapes,
    design = design,
    counts = counts,
    R_breaks = R_breaks
  )
}

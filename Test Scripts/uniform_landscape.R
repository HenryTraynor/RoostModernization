R_min <- 0.2
R_max <- 1.2

n_R_bins <- 12

R_breaks <- seq(
  R_min,
  R_max,
  length.out = n_R_bins + 1
)

R_breaks

generate_uniform_R_landscapes <- function(
    n_per_bin,
    R_breaks,
    num_roosts,
    prop_modern = 0.5,
    num_clusters_values,
    sd_values,
    max_attempts = 10000
) {
  
  n_bins <- length(R_breaks) - 1
  
  accepted <- vector("list", n_bins)
  counts <- integer(n_bins)
  
  for (b in seq_len(n_bins)) {
    accepted[[b]] <- list()
  }
  
  for (attempt in seq_len(max_attempts)) {
    
    if (all(counts >= n_per_bin)) {
      break
    }
    
    # Generate a candidate landscape
    nc <- sample(num_clusters_values, 1)
    cluster_sd <- sample(sd_values, 1)
    
    roost_info <- create_roosts(
      num_roosts = num_roosts,
      num_clusters = nc,
      sd = cluster_sd,
      prop_modern = prop_modern
    )
    
    R_value <- calculate_clark_evans_R(
      roost_info$x,
      roost_info$y
    )
    
    # Identify the target interval
    bin <- findInterval(
      R_value,
      R_breaks,
      rightmost.closed = TRUE
    )
    
    if (bin < 1 || bin > n_bins) {
      next
    }
    
    # Keep this landscape if its bin still needs samples
    if (counts[bin] < n_per_bin) {
      
      counts[bin] <- counts[bin] + 1
      
      accepted[[bin]][[counts[bin]]] <- list(
        roost_info = roost_info,
        Clark_Evans_R = R_value,
        num_clusters = nc,
        sd = cluster_sd
      )
    }
  }
  
  if (any(counts < n_per_bin)) {
    warning(
      "Some bins were not filled. Counts: ",
      paste(counts, collapse = ", "),
      ". Consider expanding the candidate parameter ranges."
    )
  }
  
  list(
    landscapes = accepted,
    counts = counts,
    R_breaks = R_breaks
  )
}

landscape_set <- generate_uniform_R_landscapes(
  n_per_bin = 10,
  R_breaks = R_breaks,
  num_roosts = 100,
  num_clusters_values = spatial_distribution$num_clusters,
  sd_values = spatial_distribution$sd
)

R_values <- unlist(
  lapply(
    landscape_set$landscapes,
    function(bin) {
      vapply(
        bin,
        function(x) x$Clark_Evans_R,
        numeric(1)
      )
    }
  )
)

hist(
  R_values,
  breaks = R_breaks,
  main = "Realized Clark-Evans R distribution",
  xlab = "Clark-Evans R"
)

table(
  cut(
    R_values,
    breaks = R_breaks,
    include.lowest = TRUE
  )
)

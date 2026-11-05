#' @importFrom stats runif rnorm dist
NULL

#' Create a set of roosts
#'
#' Creates spatially clustered roosts with a prescribed proportion of
#' modern and traditional roosts, each with different carrying capacities.
#'
#' @param num_roosts Number of roosts.
#' @param num_clusters Number of spatial clusters.
#' @param sd Standard deviation of roost locations around cluster centers.
#' @param prop_modern Proportion of roosts that are modern, between 0 and 1.
#' @param min_N_max Minimum roost carrying capacity.
#'
#' @return A data frame containing roost locations, types, carrying
#'   capacities, and population compartments.
#'
#' @export
create_roosts <- function(
    num_roosts = 25,
    num_clusters = 5,
    sd = 0.05,
    prop_modern = 0.5,
    min_N_max = 10
) {

  # Validate inputs
  stopifnot(
    num_roosts > 0,
    num_clusters > 0,
    prop_modern >= 0,
    prop_modern <= 1
  )

  # Assign roost types
  n_modern <- round(num_roosts * prop_modern)

  roost_type <- sample(
    c(
      rep("modern", n_modern),
      rep("traditional", num_roosts - n_modern)
    )
  )

  # Set mean carrying capacity by roost type
  mean_N_max <- ifelse(
    roost_type == "modern",
    100,
    20
  )

  # Cluster centers
  cluster_x <- runif(num_clusters, min = 0, max = 1)
  cluster_y <- runif(num_clusters, min = 0, max = 1)

  # Assign roosts to clusters
  cluster_id <- sample(
    seq_len(num_clusters),
    size = num_roosts,
    replace = TRUE
  )

  # Generate roost locations around cluster centers
  x <- cluster_x[cluster_id] +
    rnorm(num_roosts, mean = 0, sd = sd)

  y <- cluster_y[cluster_id] +
    rnorm(num_roosts, mean = 0, sd = sd)

  # Keep roosts within the landscape
  x <- pmin(pmax(x, 0), 1)
  y <- pmin(pmax(y, 0), 1)

  # Generate carrying capacities
  N_max <- pmax(
    round(rnorm(
      num_roosts,
      mean = mean_N_max,
      sd = mean_N_max * 0.25
    )),
    min_N_max
  )

  data.frame(
    roost_id = seq_len(num_roosts),
    roost_type = roost_type,
    x = x,
    y = y,
    N_max = N_max,
    S = 0,
    I = 0,
    R = 0,
    N = 0
  )
}

#' Initialize roost populations
#'
#' Assigns an initial population to each roost and divides individuals
#' between susceptible and infected compartments.
#'
#' @param roost_info Data frame containing roost capacities.
#' @param initial_occupancy Proportion of capacity initially occupied.
#'
#' @return Updated roost information data frame.
#'
#' @export
initialize_roost_populations <- function(
    roost_info,
    initial_occupancy = 0.25
) {

  N <- round(
    roost_info$N_max * initial_occupancy
  )

  roost_info$N <- N
  roost_info$S <- N
  roost_info$I <- 0
  roost_info$R <- 0

  roost_info
}

#' Introduce initial infections into the landscape
#'
#' Infects a specified roost or randomly selects a roost with
#' susceptible bats if no roost is specified.
#'
#' @param roost_info Initialized roost data frame.
#' @param initial_infected Number of bats to infect.
#' @param roost_id Optional index of the roost to infect. If NULL,
#'   a random roost with susceptible bats is selected.
#'
#' @return Updated roost data frame.
#' @export
seed_infection <- function(
    roost_info,
    initial_infected = 1,
    roost_id = NULL
) {

  # Identify roosts that have enough susceptible bats
  eligible_roosts <- which(roost_info$S >= initial_infected)

  if (length(eligible_roosts) == 0) {
    stop("No roost has enough susceptible bats to seed the infection.")
  }

  # Randomly select an eligible roost if none is specified
  if (is.null(roost_id)) {
    roost_id <- sample(eligible_roosts, size = 1)
  } else {

    if (!(roost_id %in% eligible_roosts)) {
      stop("The selected roost does not have enough susceptible bats.")
    }
  }

  # Move bats from susceptible to infected
  roost_info$S[roost_id] <-
    roost_info$S[roost_id] - initial_infected

  roost_info$I[roost_id] <-
    roost_info$I[roost_id] + initial_infected

  roost_info$N[roost_id] <-
    roost_info$S[roost_id] +
    roost_info$I[roost_id] +
    roost_info$R[roost_id]

  roost_info
}

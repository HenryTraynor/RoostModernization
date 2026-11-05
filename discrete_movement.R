calculate_isolation <- function(roost_info, lambda = 5) {
  
  distances <- as.matrix(
    dist(roost_info[, c("x", "y")])
  )
  
  connectivity <- exp(-lambda * distances)
  diag(connectivity) <- 0
  
  connectivity <- rowSums(connectivity)
  
  if (diff(range(connectivity)) == 0) {
    isolation <- rep(0.5, length(connectivity))
  } else {
    isolation <- 1 - (
      (connectivity - min(connectivity)) /
        diff(range(connectivity))
    )
  }
  
  isolation
}

calculate_stay_prob <- function(
    isolation,
    beta_0 = 0,
    beta_1 = 2
) {
  
  plogis(beta_0 + beta_1 * isolation)
  
}


#' Calculate roost-to-roost movement probabilities
#'
#' Creates a transition matrix based on inter-roost distance,
#' destination attractiveness, and the probability of staying.
#'
#' @param roost_info Data frame containing x, y, N, and N_max.
#' @param lambda Strength of distance decay.
#' @param alpha Strength of population-based attraction.
#' @param stay_prob Probability of remaining at the origin roost.
#' @param epsilon Small positive constant preventing zero attraction.
#'
#' @return A row-stochastic movement probability matrix.
#' @export

make_movement_matrix <- function(
    roost_info,
    lambda = 5,
    alpha = 1,
    stay_prob,
    epsilon = 0.01
) {
  
  n_roosts <- nrow(roost_info)
  
  # Allow a single stay probability or one per roost
  if (length(stay_prob) == 1) {
    stay_prob <- rep(stay_prob, n_roosts)
  }
  
  # Calculate distances between roosts
  coordinates <- as.matrix(roost_info[, c("x", "y")])
  dist_matrix <- as.matrix(stats::dist(coordinates))
  
  # Convert distances into movement weights
  connectivity_matrix <- exp(-lambda * dist_matrix)
  
  # Bats cannot move to their own roost through the movement weights;
  # staying is handled separately below.
  diag(connectivity_matrix) <- 0
  
  # Calculate attractiveness of each destination roost
  attractiveness <- (roost_info$N + epsilon)^alpha *
    pmax(1 - roost_info$N / roost_info$N_max, 0)
  
  # Weight destinations by their attractiveness
  movement_weights <- sweep(
    connectivity_matrix,
    MARGIN = 2,
    STATS = attractiveness,
    FUN = "*"
  )
  
  movement_matrix <- matrix(
    0,
    nrow = n_roosts,
    ncol = n_roosts
  )
  
  for (i in seq_len(n_roosts)) {
    weights <- movement_weights[i, ]
    total_weight <- sum(weights)
    
    if (total_weight <= 0) {
      # No available destination: remain at the origin.
      movement_matrix[i, i] <- 1
    } else {
      movement_matrix[i, ] <-
        (1 - stay_prob) * weights / total_weight
      
      movement_matrix[i, i] <- stay_prob
    }
  }
  
  dimnames(movement_matrix) <- list(
    origin = seq_len(n_roosts),
    destination = seq_len(n_roosts)
  )
  
  movement_matrix
}



#' Assign bats from one roost to destination roosts
#'
#' Samples destination roosts using a movement probability matrix,
#' while preventing roost capacities from being exceeded.
#'
#' @param n_bats Number of bats to assign.
#' @param origin Index of the origin roost.
#' @param movement_matrix Roost-to-roost movement probabilities.
#' @param N_current Current total occupancy of each destination roost.
#' @param N_max Maximum capacity of each roost.
#'
#' @return Integer vector of destination counts.
#' @export
assign_roost_movements <- function(
    n_bats,
    origin,
    movement_matrix,
    N_current,
    N_max
) {
  n_roosts <- length(N_current)
  
  assignments <- integer(n_roosts)
  
  if (n_bats == 0) {
    return(assignments)
  }
  
  capacity_remaining <- pmax(N_max - N_current, 0)
  
  for (bat in seq_len(n_bats)) {
    available <- which(capacity_remaining > 0)
    
    if (length(available) == 0) {
      stop("Insufficient roost capacity to assign all bats.")
    }
    
    probabilities <- movement_matrix[origin, available]
    probabilities <- pmax(probabilities, 0)
    
    if (sum(probabilities) <= 0) {
      # Fallback to any roost with available capacity
      probabilities <- rep(1, length(available))
    }
    
    probabilities <- probabilities / sum(probabilities)
    
    destination <- sample(
      available,
      size = 1,
      prob = probabilities
    )
    
    assignments[destination] <- assignments[destination] + 1L
    capacity_remaining[destination] <-
      capacity_remaining[destination] - 1L
  }
  
  assignments
}


#' Simulate overnight movement between roosts
#'
#' Moves susceptible, infected, and recovered bats between roosts,
#' respecting destination capacities.
#'
#' @param roost_info Data frame containing roost populations and capacities.
#' @param movement_params List containing lambda, alpha, and stay_prob.
#'
#' @return List containing updated roost populations and the movement matrix.
#' @export
run_night <- function(
    roost_info,
    movement_params = list()
) {
  lambda <- movement_params$lambda %||% 5
  alpha <- movement_params$alpha %||% 1
  stay_prob <- movement_params$stay_prob %||% 0.5
  epsilon <- movement_params$epsilon %||% 0.01
  
  required_columns <- c("x", "y", "S", "I", "R", "N", "N_max")
  
  movement_matrix <- make_movement_matrix(
    roost_info = roost_info,
    lambda = lambda,
    alpha = alpha,
    stay_prob = stay_prob,
    epsilon = epsilon
  )
  
  n_roosts <- nrow(roost_info)
  N_current <- integer(n_roosts)
  
  S_new <- integer(n_roosts)
  I_new <- integer(n_roosts)
  R_new <- integer(n_roosts)
  
  # Assign each compartment from each origin roost.
  # Update shared occupancy after every assignment.
  for (origin in seq_len(n_roosts)) {
    S_assignment <- assign_roost_movements(
      n_bats = roost_info$S[origin],
      origin = origin,
      movement_matrix = movement_matrix,
      N_current = N_current,
      N_max = roost_info$N_max
    )
    S_new <- S_new + S_assignment
    N_current <- N_current + S_assignment
    
    I_assignment <- assign_roost_movements(
      n_bats = roost_info$I[origin],
      origin = origin,
      movement_matrix = movement_matrix,
      N_current = N_current,
      N_max = roost_info$N_max
    )
    I_new <- I_new + I_assignment
    N_current <- N_current + I_assignment
    
    R_assignment <- assign_roost_movements(
      n_bats = roost_info$R[origin],
      origin = origin,
      movement_matrix = movement_matrix,
      N_current = N_current,
      N_max = roost_info$N_max
    )
    R_new <- R_new + R_assignment
    N_current <- N_current + R_assignment
  }
  
  roosts_next <- roost_info
  roosts_next$S <- S_new
  roosts_next$I <- I_new
  roosts_next$R <- R_new
  roosts_next$N <- S_new + I_new + R_new
  
  stopifnot(
    all(roosts_next$N <= roosts_next$N_max),
    sum(roosts_next$S) == sum(roost_info$S),
    sum(roosts_next$I) == sum(roost_info$I),
    sum(roosts_next$R) == sum(roost_info$R)
  )
  
  list(
    roosts = roosts_next,
    movement_matrix = movement_matrix
  )
}

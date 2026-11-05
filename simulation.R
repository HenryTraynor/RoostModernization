#' Run one complete simulation day
#'
#' Applies daytime epidemiology followed by nighttime roost movement.
#'
#' @param roost_info Data frame containing the current roost states.
#' @param beta Transmission rate.
#' @param gamma Recovery rate.
#' @param model Epidemiological model: "SIR" or "SIRS".
#' @param omega Waning immunity rate, required for SIRS.
#' @param movement_params Named list of parameters passed to run_night().
#'
#' @return A list containing updated roost states and the movement matrix.
#' @export
run_day <- function(
    roost_info,
    beta,
    gamma,
    model = c("SIR", "SIRS"),
    omega = NULL,
    movement_params = list()
) {
  model <- match.arg(model)
  
  # Run daytime epidemiology
  roost_day <- run_daytime(
    roost_info = roost_info,
    beta = beta,
    gamma = gamma,
    model = model,
    omega = omega
  )
  
  # Run nighttime movement
  night_result <- run_night(
    roost_info = roost_day,
    movement_params = movement_params
  )
  
  list(
    roosts = night_result$roosts,
    movement_matrix = night_result$movement_matrix
  )
}



#' Run a multi-day spatial epidemic simulation
#'
#' @param roost_info Initialized roost data frame. It can start without infection.
#' @param n_days Number of days to simulate.
#' @param beta Transmission rate.
#' @param gamma Recovery rate.
#' @param model Epidemiological model: "SIR" or "SIRS".
#' @param omega Waning immunity rate, required for SIRS.
#' @param movement_params Named list passed to run_night().
#' @param infection_day Day on which infection is introduced (0 to n_days).
#' @param initial_infected Number of bats initially infected.
#' @param infection_roost_id Optional row index of the roost to seed.
#'
#' @return A list containing roost_history, movement_history, and parameters.
#' @export
run_simulation <- function(
    roost_info,
    n_days = 150,
    beta,
    gamma,
    model = c("SIR", "SIRS"),
    omega = NULL,
    movement_params = list(),
    progress_every = 10,
    infection_day = 0,
    initial_infected = 1,
    infection_roost_id = NULL
) {
  model <- match.arg(model)
  
  if (model == "SIRS" && is.null(omega)) {
    stop("omega must be supplied when model = 'SIRS'.")
  }
  
  # Save day 0, before any infection if infection_day > 0
  roost_history <- vector("list", n_days + 1L)
  movement_history <- vector("list", n_days)
  roost_history[[1]] <- roost_info
  
  for (day in seq_len(n_days)) {
    
    # Seed infection at the beginning of the chosen day.
    # Day 0 is handled separately below.
    if (day == infection_day) {
      roost_history[[day]] <- seed_infection(
        roost_history[[day]],
        initial_infected = initial_infected,
        roost_id = infection_roost_id
      )
    }
    
    result <- run_day(
      roost_info = roost_history[[day]],
      beta = beta,
      gamma = gamma,
      model = model,
      omega = omega,
      movement_params = movement_params
    )
    
    roost_history[[day + 1L]] <- result$roosts
    movement_history[[day]] <- result$movement_matrix
  }
  
  # If infection_day = 0, seed before the first simulated day.
  if (infection_day == 0) {
    roost_history[[1]] <- seed_infection(
      roost_info,
      initial_infected = initial_infected,
      roost_id = infection_roost_id
    )
  }
  
  names(roost_history) <- paste0("day_", 0:n_days)
  names(movement_history) <- paste0("night_", seq_len(n_days))
  
  list(
    roost_history = roost_history,
    movement_history = movement_history,
    parameters = list(
      n_days = n_days,
      beta = beta,
      gamma = gamma,
      model = model,
      omega = omega,
      movement_params = movement_params,
      infection_day = infection_day,
      initial_infected = initial_infected,
      infection_roost_id = infection_roost_id
    )
  )
}

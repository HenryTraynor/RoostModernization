#' @importFrom stats rbinom
NULL

#' Perform one SIR epidemiological time step
#'
#' Updates susceptible, infected, and recovered individuals in a roost using binomial transitions for infection and recovery.
#'
#' @param y Numeric vector containing S, I, and R.
#' @param beta Transmission rate.
#' @param gamma Recovery rate.
#'
#' @return Numeric vector containing updated S, I, and R.
#'
#' @export
SIR_step <- function(y, beta, gamma) {

  S <- unname(y[1])
  I <- unname(y[2])
  R <- unname(y[3])

  N <- S + I + R

  # No infection possible if there are no individuals
  if (N <= 0) {
    return(c(S = S, I = I, R = R))
  }

  # Probability of infection during this time step
  p_infection <- 1 - exp(-beta * I / N)

  # Probability of recovery during this time step
  p_recovery <- 1 - exp(-gamma)

  # Number of susceptible individuals becoming infected
  new_infections <- rbinom(
    1,
    size = S,
    prob = p_infection
  )

  # Number of infected individuals recovering
  new_recoveries <- rbinom(
    1,
    size = I,
    prob = p_recovery
  )

  # Update compartments
  S_new <- S - new_infections
  I_new <- I + new_infections - new_recoveries
  R_new <- R + new_recoveries

  c(
    S = S_new,
    I = I_new,
    R = R_new
  )
}


#' Perform one SIRS epidemiological time step
#'
#' Updates susceptible, infected, and recovered individuals in a roost,
#' including loss of immunity from the recovered compartment.
#'
#' @param y Numeric vector containing S, I, and R.
#' @param beta Transmission rate.
#' @param gamma Recovery rate.
#' @param omega Rate of loss of immunity.
#'
#' @return Numeric vector containing updated S, I, and R.
#'
#' @export
SIRS_step <- function(y, beta, gamma, omega) {

  S <- unname(y[1])
  I <- unname(y[2])
  R <- unname(y[3])

  N <- S + I + R

  if (N <= 0) {
    return(c(S = S, I = I, R = R))
  }

  # Transition probabilities
  p_infection <- 1 - exp(-beta * I / N)
  p_recovery <- 1 - exp(-gamma)
  p_loss_immunity <- 1 - exp(-omega)

  # Draw transitions
  new_infections <- rbinom(
    1,
    size = S,
    prob = p_infection
  )

  new_recoveries <- rbinom(
    1,
    size = I,
    prob = p_recovery
  )

  loss_immunity <- rbinom(
    1,
    size = R,
    prob = p_loss_immunity
  )

  # Update compartments
  S_new <- S - new_infections + loss_immunity

  I_new <- I +
    new_infections -
    new_recoveries

  R_new <- R +
    new_recoveries -
    loss_immunity

  c(
    S = S_new,
    I = I_new,
    R = R_new
  )
}


#' Run daytime epidemiology across roosts
#'
#' Applies an SIR or SIRS model independently to each roost.
#'
#' @param roost_info Data frame containing S, I, and R populations.
#' @param beta Transmission rate.
#' @param gamma Recovery rate.
#' @param model Epidemiological model. Either `"SIR"` or `"SIRS"`.
#' @param omega Rate of loss of immunity. Required for the SIRS model.
#'
#' @return Updated roost information data frame.
#'
#' @export
run_daytime <- function(
    roost_info,
    beta,
    gamma,
    model = "SIR",
    omega = NULL
) {

  model <- match.arg(model, c("SIR", "SIRS"))

  if (model == "SIRS" && is.null(omega)) {
    stop("omega must be provided when model = 'SIRS'.")
  }

  roosts <- roost_info

  for (i in seq_len(nrow(roosts))) {

    y <- c(
      S = roosts$S[i],
      I = roosts$I[i],
      R = roosts$R[i]
    )

    if (model == "SIR") {

      y_new <- SIR_step(
        y = y,
        beta = beta,
        gamma = gamma
      )

    } else {

      y_new <- SIRS_step(
        y = y,
        beta = beta,
        gamma = gamma,
        omega = omega
      )
    }

    roosts$S[i] <- y_new["S"]
    roosts$I[i] <- y_new["I"]
    roosts$R[i] <- y_new["R"]

    roosts$N[i] <-
      roosts$S[i] +
      roosts$I[i] +
      roosts$R[i]
  }

  roosts
}

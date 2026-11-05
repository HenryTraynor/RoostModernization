
#' Run daytime epidemiological dynamics across roosts
#'
#' Applies an SIR or SIRS update independently to each roost.
#'
#' @param roost_info Data frame containing S, I, R, and N.
#' @param beta Transmission rate.
#' @param gamma Recovery rate.
#' @param model Epidemiological model: "SIR" or "SIRS".
#' @param omega Rate of waning immunity, required for SIRS.
#'
#' @return Updated roost data frame.
#' @export
run_daytime <- function(
    roost_info,
    beta,
    gamma,
    model = c("SIR", "SIRS"),
    omega = NULL
) {
  model <- match.arg(model)
  
  if (model == "SIRS" && is.null(omega)) {
    stop("omega must be supplied when model = 'SIRS'.")
  }
  
  for (i in seq_len(nrow(roost_info))) {
    y <- c(
      S = roost_info$S[i],
      I = roost_info$I[i],
      R = roost_info$R[i]
    )
    
    if (model == "SIR") {
      y_new <- SIR_step(y, beta = beta, gamma = gamma)
    } else {
      y_new <- SIRS_step(
        y,
        beta = beta,
        gamma = gamma,
        omega = omega
      )
    }
    
    roost_info$S[i] <- y_new["S"]
    roost_info$I[i] <- y_new["I"]
    roost_info$R[i] <- y_new["R"]
    roost_info$N[i] <- sum(y_new[c("S", "I", "R")])
  }
  
  roost_info
}

source("epidemiology.R")
source("initialize_roosts.R")

set.seed(123)

#------------
# STEP FUNCTIONS
#------------
# Helper to check compartment counts
check_counts <- function(y) {
  stopifnot(
    all(c("S", "I", "R") %in% names(y)),
    all(is.finite(y[c("S", "I", "R")])),
    all(y[c("S", "I", "R")] >= 0),
    all(y[c("S", "I", "R")] == floor(y[c("S", "I", "R")]))
  )
}

# Initial population
y0 <- c(S = 95, I = 5, R = 0)
N0 <- sum(y0)

# Test SIR
sir_result <- SIR_step(y0, beta = 0.4, gamma = 0.2)

check_counts(sir_result)
stopifnot(sum(sir_result[c("S", "I", "R")]) == N0)

print(sir_result)

# Test SIRS
sirs_result <- SIRS_step(
  y0,
  beta = 0.4,
  gamma = 0.2,
  omega = 0.05
)

check_counts(sirs_result)
stopifnot(sum(sirs_result[c("S", "I", "R")]) == N0)

print(sirs_result)


#---------------
# INF. SEEDING
#---------------

roosts <- create_roosts(
  num_roosts = 10,
  num_clusters = 3,
  prop_modern = 0.5
)

roosts <- initialize_roost_populations(
  roosts,
  initial_occupancy = 0.25
)

N_before <- sum(roosts$N)
I_before <- sum(roosts$I)

# Seed 3 infected bats in roost 1
roosts_seeded <- seed_infection(
  roosts,
  initial_infected = 3,
  roost_id = 1
)

stopifnot(
  sum(roosts_seeded$N) == N_before,
  sum(roosts_seeded$I) == I_before + 3,
  roosts_seeded$S[1] == roosts$S[1] - 3,
  roosts_seeded$I[1] == roosts$I[1] + 3,
  all(roosts_seeded$N ==
        roosts_seeded$S + roosts_seeded$I + roosts_seeded$R)
)

message("seed_infection() passed!")


set.seed(42)

roosts_random <- seed_infection(
  roosts,
  initial_infected = 1
)

which(roosts_random$I > roosts$I)

#--------------------
# VISUALIZATION
#------------------

simulate_epidemic <- function(
    model = c("SIR", "SIRS"),
    n_days = 150,
    initial_population = 100,
    initial_infected = 5,
    beta = 0.4,
    gamma = 0.2,
    omega = 0.05
) {
  model <- match.arg(model)
  
  y <- c(
    S = initial_population - initial_infected,
    I = initial_infected,
    R = 0
  )
  
  history <- data.frame(
    day = 0,
    S = y["S"],
    I = y["I"],
    R = y["R"],
    model = model
  )
  
  for (day in seq_len(n_days)) {
    if (model == "SIR") {
      y <- SIR_step(y, beta = beta, gamma = gamma)
    } else {
      y <- SIRS_step(
        y,
        beta = beta,
        gamma = gamma,
        omega = omega
      )
    }
    
    history <- rbind(
      history,
      data.frame(
        day = day,
        S = y["S"],
        I = y["I"],
        R = y["R"],
        model = model
      )
    )
  }
  
  history
}

set.seed(123)

sir_history <- simulate_epidemic(model = "SIR",
                                 n_days = 300)
sirs_history <- simulate_epidemic(model = "SIRS",
                                  n_days = 300)

epidemic_history <- rbind(sir_history, sirs_history)

library(ggplot2)
library(tidyr)

epidemic_long <- pivot_longer(
  epidemic_history,
  cols = c(S, I, R),
  names_to = "compartment",
  values_to = "bats"
)

ggplot(
  epidemic_long,
  aes(x = day, y = bats, color = compartment)
) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~model, ncol = 1) +
  labs(
    title = "SIR vs. SIRS epidemic dynamics",
    x = "Day",
    y = "Number of bats",
    color = "Compartment"
  ) +
  theme_minimal(base_size = 12)

#------------------
# Many sims
#-----------------

run_replicates <- function(n_replicates = 100, n_days = 150) {
  results <- vector("list", n_replicates * 2)
  index <- 1
  
  for (model in c("SIR", "SIRS")) {
    for (rep in seq_len(n_replicates)) {
      results[[index]] <- simulate_epidemic(
        model = model,
        n_days = n_days,
        initial_population = 100,
        initial_infected = 5,
        beta = 0.4,
        gamma = 0.2,
        omega = 0.05
      ) |>
        transform(replicate = rep)
      
      index <- index + 1
    }
  }
  
  do.call(rbind, results)
}

set.seed(123)

replicate_results <- run_replicates(n_replicates = 100)

infected_summary <- aggregate(
  I ~ model + day,
  data = replicate_results,
  FUN = function(x) {
    c(
      mean = mean(x),
      lower = unname(quantile(x, 0.025)),
      upper = unname(quantile(x, 0.975))
    )
  }
)

infected_summary <- data.frame(
  model = infected_summary$model,
  day = infected_summary$day,
  infected_summary$I
)

# Check the result
head(infected_summary)

ggplot(
  infected_summary,
  aes(x = day, y = mean, color = model, fill = model)
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.2,
    color = NA
  ) +
  geom_line(linewidth = 1) +
  labs(
    title = "Mean epidemic trajectory across simulations",
    subtitle = "Shaded regions show 95% simulation intervals",
    x = "Day",
    y = "Infected bats",
    color = "Model",
    fill = "Model"
  ) +
  theme_minimal(base_size = 12)


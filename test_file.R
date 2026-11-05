landscape_set <- generate_landscapes(
  n_per_bin = 5,
  R_breaks = c(0, 0.4, 0.8, 1.2),
  prop_modern_values = seq(0, 1, by = 0.2),
  num_roosts = 30,
  cluster_values = c(3, 5, 8),
  sd_values = c(0.03, 0.05, 0.1)
)

# Number of landscape scenarios
design <- landscape_set$design
n_sims <- nrow(design)

# Storage for simulation results
simulation_results <- vector("list", n_sims)

# Reproducibility
set.seed(123)

for (i in seq_len(n_sims)) {
  
  message("Running simulation ", i, " of ", n_sims)
  
  # Extract landscape
  roost_info <- design$roost_info[[i]]
  
  # Initialize bat populations
  roost_info <- initialize_roost_populations(
    roost_info,
    initial_occupancy = 0.25
  )
  
  # Run epidemic simulation
  simulation_results[[i]] <- run_simulation(
    roost_info = roost_info,
    n_days = 150,
    beta = 0.4,
    gamma = 0.2,
    model = "SIR",
    infection_day = 30,
    initial_infected = 5,
    movement_params = list(
      lambda = 5,
      alpha = 1,
      stay_prob = 0.5,
      epsilon = 0.01
    )
  )
  
  # Attach landscape information
  simulation_results[[i]]$landscape_info <- list(
    simulation_id = i,
    landscape_id = design$landscape_id[i],
    prop_modern = design$prop_modern[i],
    Clark_Evans_R = design$Clark_Evans_R[i],
    num_clusters = design$num_clusters[i],
    sd = design$sd[i],
    R_bin = design$R_bin[i]
  )
}

names(simulation_results) <- paste0(
  "simulation_",
  seq_len(n_sims)
)

library(RoostEpiPkg)

set.seed(123)

roost_info <- create_roosts(
  num_roosts = 25,
  num_clusters = 5,
  prop_modern = 0.5
)

roost_info <- initialize_roost_populations(
  roost_info,
  initial_occupancy = 0.25,
  initial_infected = 0.25
)

spatial_params <- create_spatial_params(nx = 100, ny = 100)

movement_params <- setup_movement(
  roost_info = roost_info,
  spatial_params = spatial_params
)

epidemic_params <- create_epidemic_params(
  beta = 0.5,
  gamma = 0.03,
  model = "SIR"
)

sim <- run_simulation(
  roost_info = roost_info,
  spatial_params = spatial_params,
  movement_params = movement_params,
  epidemic_params = epidemic_params,
  n_days = 5
)

str(sim, max.level = 2)
names(sim[[1]])
head(sim[[1]])

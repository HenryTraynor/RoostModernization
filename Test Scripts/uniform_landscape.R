R_min <- 0.2
R_max <- 1.2

n_R_bins <- 12

R_breaks <- seq(
  R_min,
  R_max,
  length.out = n_R_bins + 1
)

R_breaks

landscape_set <- generate_landscapes(
  n_per_bin = 10,
  R_breaks = R_breaks,
  prop_modern_values = seq(0, 1, by = 0.2),
  num_roosts = 100,
  cluster_values = cluster_seq,
  sd_values = sd_seq,
  modern_N_max = 100,
  traditional_N_max = 20
)

simulation_grid <- landscape_set$design

nrow(simulation_grid)
# Expected: 450

design <- landscape_set$design

r25 <- design$roost_info[[
  which(
    design$landscape_index == 1 &
      design$prop_modern == 0.2
  )
]]

r50 <- design$roost_info[[
  which(
    design$landscape_index == 1 &
      design$prop_modern == 0.4
  )
]]

# Coordinates remain fixed
stopifnot(
  identical(r25$x, r50$x),
  identical(r25$y, r50$y)
)

# Any roost that stays the same type retains its capacity
same_type <- r25$roost_type == r50$roost_type

stopifnot(
  all(r25$N_max[same_type] == r50$N_max[same_type])
)




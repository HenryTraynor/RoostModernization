set.seed(123)

landscape_set <- generate_landscapes(
  n_per_bin = 5,
  R_breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2),
  prop_modern_values = c(0, 0.5, 1),
  num_roosts = 30,
  cluster_values = c(3, 5, 8),
  sd_values = c(0.03, 0.05, 0.1),
  max_attempts = 10000
)

landscape_set$counts
head(landscape_set$design[, -which(names(landscape_set$design) == "roost_info")])

# Check that each scenario has the same coordinates as its base landscape
landscape_set$design$roost_info[[1]][, c("x", "y", "roost_type", "N_max")]

design <- landscape_set$design

# Select one base landscape
landscape_id <- unique(design$landscape_id)[1]

# Get all modernization scenarios for that landscape
scenario_rows <- which(design$landscape_id == landscape_id)

# Extract roost dataframes
scenario_roosts <- design$roost_info[scenario_rows]

# Compare potential modern capacities
modern_capacities <- lapply(
  scenario_roosts,
  function(x) {
    # Reconstructing potential capacities requires
    # the original saved base landscape.
    x$N_max
  }
)

# Display the actual assigned capacities
for (i in seq_along(scenario_rows)) {
  cat(
    "\nProportion modern:",
    design$prop_modern[scenario_rows[i]],
    "\n"
  )
  
  print(data.frame(
    roost_id = scenario_roosts[[i]]$roost_id,
    type = scenario_roosts[[i]]$roost_type,
    N_max = scenario_roosts[[i]]$N_max
  ))
}


base_landscape <- landscape_set$landscapes[[landscape_id]]

modernization_order <- base_landscape$modernization_order

# Check each scenario
for (i in seq_along(scenario_rows)) {
  
  scenario <- scenario_roosts[[i]]
  
  expected_n_modern <- round(
    design$prop_modern[scenario_rows[i]] *
      nrow(scenario)
  )
  
  expected_modern <- head(
    modernization_order,
    expected_n_modern
  )
  
  actual_modern <- which(scenario$roost_type == "modern")
  
  stopifnot(setequal(expected_modern, actual_modern))
}

message("Modernization order is consistent!")

modern_N_max <- 100
traditional_N_max <- 20

base_roosts <- base_landscape$base_roosts

summary(base_roosts$N_max_modern)
summary(base_roosts$N_max_traditional)

data.frame(
  building_type = c("Modern", "Traditional"),
  target_mean = c(modern_N_max, traditional_N_max),
  actual_mean = c(
    mean(base_roosts$N_max_modern),
    mean(base_roosts$N_max_traditional)
  ),
  target_sd = c(
    0.25 * modern_N_max,
    0.25 * traditional_N_max
  ),
  actual_sd = c(
    sd(base_roosts$N_max_modern),
    sd(base_roosts$N_max_traditional)
  )
)


library(ggplot2)
library(tidyr)

capacity_data <- data.frame(
  Traditional = base_roosts$N_max_traditional,
  Modern = base_roosts$N_max_modern
) |>
  pivot_longer(
    cols = everything(),
    names_to = "building_type",
    values_to = "N_max"
  )

ggplot(
  capacity_data,
  aes(x = N_max, fill = building_type)
) +
  geom_histogram(
    bins = 20,
    alpha = 0.65,
    position = "identity"
  ) +
  facet_wrap(~building_type, scales = "free") +
  labs(
    title = "Potential roost capacity distributions",
    x = "Maximum roost capacity",
    y = "Number of roosts"
  ) +
  theme_minimal()

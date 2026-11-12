source("discrete_movement.R")

set.seed(123)

movement_params = list(
  lambda = 5,
  alpha = 5,
  stay_prob = 0.75,
  epsilon = 0.01
)

roost_info <- data.frame(
  x = c(0.1, 0.3, 0.5, 0.7, 0.9),
  y = c(0.2, 0.8, 0.4, 0.7, 0.3),
  S = c(15, 10, 20, 8, 12),
  I = c(3, 2, 5, 1, 4),
  R = c(2, 3, 5, 1, 4),
  N_max = c(50, 10, 10, 50, 50)
)

roost_info$N <- with(roost_info, S + I + R)

print(roost_info)

movement_matrix <- make_movement_matrix(
  roost_info = roost_info,
  lambda = movement_params$lambda,
  alpha = movement_params$alpha,
  stay_prob = movement_params$stay_prob,
  epsilon = movement_params$epsilon
)

# Simulation settings
n_days <- 365

# Store roost populations for every day
roost_history <- vector("list", n_days + 1)

# Save initial conditions (day 0)
roost_history[[1]] <- roost_info

# Epidemiological parameters
beta <- 0.4
gamma <- 0.2
model <- "SIRS"
omega <- 0.05  # Set a value if using model = "SIRS"

# Movement parameters
movement_params <- list(
  lambda = 5,
  alpha = 1,
  stay_prob = 0.5,
  epsilon = 0.01
)

# Run simulation
for (day in seq_len(n_days)) {
  
  # 1. Daytime epidemiology
  roost_info <- run_daytime(
    roost_info = roost_info,
    beta = beta,
    gamma = gamma,
    model = model,
    omega = omega
  )
  
  # 2. Overnight movement
  night_output <- run_night(
    roost_info = roost_info,
    movement_params = movement_params
  )
  
  # 3. Update roost populations after movement
  roost_info <- night_output$roosts
  
  # 4. Save the state for this day
  roost_history[[day + 1]] <- roost_info
  
  # 5. Progress message
  if (day %% 10 == 0 || day == n_days) {
    message("Completed day ", day, " of ", n_days)
  }
}

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(zoo)

# Combine all simulation days into one dataframe
occupancy_history <- bind_rows(
  lapply(seq_along(roost_history), function(i) {
    dat <- roost_history[[i]]
    dat$Roost <- seq_len(nrow(dat))
    dat$Day <- i - 1
    dat
  })
)

occupancy_history <- occupancy_history %>%
  group_by(Roost) %>%
  arrange(Day) %>%
  mutate(
    N_smooth = zoo::rollmean(
      N,
      k = 5,
      fill = NA,
      align = "center"
    )
  ) %>%
  ungroup()

roost_infection_totals <- occupancy_history %>%
  group_by(Roost) %>%
  summarise(
    Total_Infected = sum(I),
    .groups = "drop"
  )

# Landscape-wide epidemic totals
epidemic_history <- occupancy_history %>%
  group_by(Day) %>%
  summarise(
    S = sum(S),
    I = sum(I),
    R = sum(R),
    N = sum(N),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(S, I, R),
    names_to = "Stage",
    values_to = "Bats"
  )

ui <- fluidPage(
  titlePanel("Bat Roost Epidemic Simulator"),
  
  sliderInput(
    "day",
    "Simulation day:",
    min = 0,
    max = length(roost_history) - 1,
    value = 0,
    step = 1,
    animate = animationOptions(
      interval = 1000,
      loop = FALSE
    )
  ),
  
  h3("Landscape-wide epidemic"),
  plotOutput("epidemic_plot", height = "300px"),
  
  fluidRow(
    column(
      width = 6,
      h3("Roost occupancy through time"),
      plotOutput("occupancy_plot", height = "400px")
    ),
    column(
      width = 6,
      h3("Spatial distribution of infection"),
      plotOutput("map_plot", height = "400px")
    )
  )
)

server <- function(input, output, session) {
  
  # Roost states for the selected day
  current_roosts <- reactive({
    roost_history[[input$day + 1]]
  })
  
  # Landscape-wide epidemic curve
  output$epidemic_plot <- renderPlot({
    
    ggplot(
      epidemic_history %>% filter(Day <= input$day),
      aes(x = Day, y = Bats, color = Stage)
    ) +
      geom_line(linewidth = 1) +
      geom_vline(
        xintercept = input$day,
        linetype = "dashed"
      ) +
      scale_color_manual(
        values = c(
          S = "#2878B5",
          I = "#D62728",
          R = "#2CA02C"
        ),
        labels = c(
          S = "Susceptible",
          I = "Infected",
          R = "Recovered"
        )
      ) +
      scale_y_continuous(labels = label_comma()) +
      labs(
        x = "Day",
        y = "Number of bats",
        color = "Epidemic stage"
      ) +
      theme_minimal(base_size = 12)
  })
  
  # Occupancy of each roost over time
  output$occupancy_plot <- renderPlot({
    
    ggplot(
      occupancy_history %>% filter(Day <= input$day),
      aes(x = Day, group = Roost, color = factor(Roost))
    ) +
      geom_line(
        aes(y = N),
        alpha = 0.12,
        linewidth = 0.4
      ) +
      geom_line(
        aes(y = N_smooth),
        linewidth = 1
      ) +
      geom_vline(
        xintercept = input$day,
        linetype = "dashed"
      ) +
      labs(
        x = "Day",
        y = "Roost occupancy",
        color = "Roost"
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none")
  })
  
  # Spatial map for the selected day
  output$map_plot <- renderPlot({
    
    roosts <- current_roosts()
    
    # Recalculate potential connectivity at the selected day's occupancy
    movement_matrix <- make_movement_matrix(
      roost_info = roosts,
      lambda = movement_params$lambda,
      alpha = movement_params$alpha,
      stay_prob = movement_params$stay_prob,
      epsilon = movement_params$epsilon
    )
    
    # Create an edge dataframe from movement probabilities
    movement_df <- as.data.frame(as.table(movement_matrix))
    names(movement_df) <- c(
      "Origin", "Destination", "Probability"
    )
    
    movement_df$Origin <- as.integer(
      as.character(movement_df$Origin)
    )
    movement_df$Destination <- as.integer(
      as.character(movement_df$Destination)
    )
    
    # Show connections with meaningful movement probability
    movement_df <- movement_df %>%
      filter(
        Origin != Destination,
        Probability >= 0.02
      ) %>%
      mutate(
        x = roosts$x[Origin],
        y = roosts$y[Origin],
        xend = roosts$x[Destination],
        yend = roosts$y[Destination]
      )
    
    # Add cumulative infections across the entire simulation
    roosts <- roosts %>%
      mutate(
        Roost = seq_len(nrow(roosts))
      ) %>%
      left_join(roost_infection_totals, by = "Roost")
    
    ggplot() +
      geom_segment(
        data = movement_df,
        aes(
          x = x, y = y,
          xend = xend, yend = yend,
          linewidth = Probability
        ),
        color = "grey55",
        alpha = 0.4,
        lineend = "round"
      ) +
      geom_point(
        data = roosts,
        aes(
          x = x,
          y = y,
          size = N,
          color = Total_Infected
        ),
        alpha = 0.9
      ) +
      geom_text(
        data = roosts,
        aes(x = x, y = y, label = round(Total_Infected)),
        nudge_y = 0.025,
        size = 3
      ) +
      scale_color_gradient(
        low = "#2C7BB6",
        high = "#D7191C",
        name = "Cumulative\ninfected bat-days"
      ) +
      scale_size_continuous(
        range = c(2, 10),
        name = "Occupancy"
      ) +
      scale_linewidth_continuous(
        range = c(0.2, 1.5),
        guide = "none"
      ) +
      coord_equal() +
      labs(
        x = "X coordinate",
        y = "Y coordinate",
        title = paste("Day", input$day)
      ) +
      theme_minimal(base_size = 12)
  })
}

shinyApp(ui = ui, server = server)


#-----------------
# STATIC VIS
#-------------------


library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(zoo)
library(patchwork)

# --------------------------------------------------
# 1. Landscape-wide epidemic
# --------------------------------------------------

epidemic_plot <- ggplot(
  epidemic_history,
  aes(x = Day, y = Bats, color = Stage)
) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(
    values = c(
      S = "#2878B5",
      I = "#D62728",
      R = "#2CA02C"
    ),
    labels = c(
      S = "Susceptible",
      I = "Infected",
      R = "Recovered"
    )
  ) +
  scale_y_continuous(labels = label_comma()) +
  labs(
    title = "Landscape-wide epidemic",
    x = "Simulation day",
    y = "Number of bats",
    color = "Epidemic stage"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )


# --------------------------------------------------
# 2. Roost occupancy through time
# --------------------------------------------------

occupancy_plot <- ggplot(
  occupancy_history,
  aes(
    x = Day,
    y = N,
    group = Roost,
    color = factor(Roost)
  )
) +
  geom_line(
    aes(y = N),
    alpha = 0.15,
    linewidth = 0.4
  ) +
  geom_line(
    aes(y = N_smooth),
    linewidth = 0.9
  ) +
  labs(
    title = "Roost occupancy through time",
    x = "Simulation day",
    y = "Roost occupancy"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )


# --------------------------------------------------
# 3. Spatial distribution of infection
# --------------------------------------------------

# Use final-day occupancy for point size
final_roosts <- roost_history[[length(roost_history)]]

final_roosts$Roost <- seq_len(nrow(final_roosts))

# Add cumulative infection burden
final_roosts <- final_roosts %>%
  left_join(
    roost_infection_totals,
    by = "Roost"
  )

# Recalculate movement probabilities at final occupancy
final_movement_matrix <- make_movement_matrix(
  roost_info = final_roosts,
  lambda = movement_params$lambda,
  alpha = movement_params$alpha,
  stay_prob = movement_params$stay_prob,
  epsilon = movement_params$epsilon
)

# Convert movement matrix to edge dataframe
movement_df <- as.data.frame(
  as.table(final_movement_matrix)
)

names(movement_df) <- c(
  "Origin",
  "Destination",
  "Probability"
)

movement_df <- movement_df %>%
  mutate(
    Origin = as.integer(as.character(Origin)),
    Destination = as.integer(as.character(Destination))
  ) %>%
  filter(
    Origin != Destination,
    Probability >= 0.02
  ) %>%
  mutate(
    x = final_roosts$x[Origin],
    y = final_roosts$y[Origin],
    xend = final_roosts$x[Destination],
    yend = final_roosts$y[Destination]
  )

# Spatial plot
map_plot <- ggplot() +
  
  geom_segment(
    data = movement_df,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      linewidth = Probability
    ),
    color = "grey55",
    alpha = 0.4,
    lineend = "round"
  ) +
  
  geom_point(
    data = final_roosts,
    aes(
      x = x,
      y = y,
      size = N,
      color = Total_Infected
    ),
    alpha = 0.9
  ) +
  
  geom_text(
    data = final_roosts,
    aes(
      x = x,
      y = y,
      label = round(Total_Infected)
    ),
    nudge_y = 0.025,
    size = 3
  ) +
  
  scale_color_gradient(
    low = "#2C7BB6",
    high = "#D7191C",
    name = "Cumulative\ninfected bat-days"
  ) +
  
  scale_size_continuous(
    range = c(3, 12),
    name = "Final occupancy"
  ) +
  
  scale_linewidth_continuous(
    range = c(0.2, 1.5),
    guide = "none"
  ) +
  
  coord_equal() +
  
  labs(
    title = "Spatial distribution of infection",
    x = "X coordinate",
    y = "Y coordinate"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )


# --------------------------------------------------
# 4. Combine visualizations
# --------------------------------------------------
library(patchwork)


combined_plot <- (
  epidemic_plot /
    occupancy_plot /
    map_plot
) +
  plot_layout(
    heights = c(1, 1.2, 1.2),
    guides = "collect"
  ) +
  plot_annotation(
    title = "Bat Roost Epidemic Simulation",
    subtitle = paste(
      "SIRS model |",
      nrow(roost_info), "roosts |",
      n_days, "simulation days"
    ),
    theme = theme(
      plot.title = element_text(
        size = 18,
        face = "bold"
      ),
      plot.subtitle = element_text(size = 11)
    )
  )

print(combined_plot)

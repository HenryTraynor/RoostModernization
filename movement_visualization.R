# Load your package functions
# source("R/distance_decay.R")
# source("R/movement.R")

set.seed(123)

movement_params = list(
  lambda = 5,
  alpha = 1,
  stay_prob = 0.75,
  epsilon = 0.01
)

roost_info <- data.frame(
  x = c(0.1, 0.3, 0.5, 0.7, 0.9),
  y = c(0.2, 0.8, 0.4, 0.7, 0.3),
  S = c(15, 10, 20, 8, 12),
  I = c(3, 2, 5, 1, 4),
  R = c(2, 3, 5, 1, 4),
  N_max = c(30, 25, 35, 20, 25)
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
n_days <- 150

# Store roost populations for every day
roost_history <- vector("list", n_days + 1)

# Save initial conditions (day 0)
roost_history[[1]] <- roost_info

# Epidemiological parameters
beta <- 0.4
gamma <- 0.2
model <- "SIR"
omega <- NULL  # Set a value if using model = "SIRS"

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

# Combine all simulation days into one dataframe
occupancy_history <- bind_rows(
  lapply(seq_along(roost_history), function(i) {
    dat <- roost_history[[i]]
    dat$Roost <- seq_len(nrow(dat))
    dat$Day <- i - 1
    dat
  })
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
      aes(x = Day, y = N, color = factor(Roost), group = Roost)
    ) +
      geom_line(alpha = 0.75) +
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
    
    # Infection prevalence within each roost
    roosts <- roosts %>%
      mutate(
        Infection = ifelse(N > 0, I / N, 0)
      )
    
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
          color = Infection
        ),
        alpha = 0.9
      ) +
      scale_color_gradient(
        low = "#2C7BB6",
        high = "#D7191C",
        limits = c(0, 1),
        labels = label_percent(),
        name = "Infected"
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

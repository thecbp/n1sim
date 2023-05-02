library(shiny)
library(tidyverse)

# Define server logic for Shiny app
server <- function(input, output) {
  
  # Storing data for plotting on dashboard
  simulated_data <- reactiveVal()
  
  ##############################################################################
  
  # REACTIVE UI
  
  ##############################################################################
  
  output$Ejs = renderUI({
    
    lapply(1:input$n_trts, function(i) {
      fluidRow(
        column(12, numericInput(paste0("ej_", i), paste0("Treatment effect ", i, ":"), value = 0))
      )
    })
    
  })
  
  output$Gammas = renderUI({
    
    lapply(1:input$n_trts, function(i) {
      fluidRow(
        column(12, numericInput(paste0("gamma_", i), paste0("Gamma ", i, "  (washout parameter):"), value = 0.2))
      )
    })
    
  })
  
  output$Taus= renderUI({
    
    lapply(1:input$n_trts, function(i) {
      fluidRow(
        column(12, numericInput(paste0("tau_", i), paste0("Tau ", i, " (run-in parameter):"), value = 0.3))
      )
    })
    
  })
  
  ##############################################################################
  
  # DERIVED VALUES
  
  ##############################################################################
  
  treatment_effects = reactive({
    
    n_trts = input$n_trts
    ejs = numeric(n_trts)
    
    for (i in 1:n_trts) {
      label = paste0("ej_", i)
      ejs[i] = input[[label]]
    }
    
    ejs
    
  })
  
  gammas = reactive({
    
    n_trts = input$n_trts
    gammas = numeric(n_trts)
    
    for (i in 1:n_trts) {
      label = paste0("gamma_", i)
      gammas[i] = input[[label]]
    }
    
    gammas
    
  })
  
  taus = reactive({
    
    n_trts = input$n_trts
    taus = numeric(n_trts)
    
    for (i in 1:n_trts) {
      label = paste0("tau_", i)
      taus[i] = input[[label]]
    }
    
    taus
    
  })
  
  
  ##############################################################################
  
  # EVENT OBSERVERS
  
  ##############################################################################
  
  observeEvent(input$simulate, {
    
    data = n1sim::simulate(n_trts = input$n_trts, 
                           n_blocks = input$n_blocks, 
                           period_length = input$period_length, 
                           washout_length = input$washout_length,
                           sampling_timestep = input$sampling_timestep, 
                           sd_b = input$sd_b, 
                           sd_p = input$sd_p, 
                           sd_o = input$sd_o,
                           Ejs = treatment_effects(), 
                           gammas = gammas(), 
                           taus = taus(), 
                           eta = input$eta, 
                           phi = input$phi)
    
    
    simulated_data(data)
    
  })
  
  
  output$test = renderDataTable({
    
    data = simulated_data()

    if (is.null(data)) {
      tibble::tibble(x = 1:3, y = 4:6)
    } else {
      data[["data"]]
    }
    
  })
  
  output$treatmentPlot = renderPlot({
    
    data = simulated_data()
    
    washout_indicator = input$n_trts + 1
    
    # Get ranges for when treatment or washout is active
    rectangles = data[["data"]] %>% 
      group_by(period) %>% 
      summarize(
        xmin = min(time),
        xmax = max(time) + 1,
        treatment = unique(treatment)
      ) %>% 
      mutate(
        Name = case_when(
          treatment == washout_indicator ~ "Washout",
          TRUE ~ paste0("X", treatment)
        ) %>% as.factor
      )
    
    effects = data[["data"]] %>% 
      select(time, matches("X[0-9]+")) %>% 
      pivot_longer(
        matches("X[0-9]$"),
        names_to = "Treatment",
        values_to = "effect"
      ) 
    
    p = ggplot() + 
      geom_line(data = effects, aes(x = time, y = effect, color = Treatment), show.legend = FALSE) +
      geom_rect(data = rectangles, aes(xmin = xmin, xmax = xmax, 
                                       ymin = -Inf, ymax = Inf, 
                                       fill = Name), 
                color = NA, alpha = 0.1) + 
      theme_minimal() +
      labs(
        x = "Observation number",
        y = "True treatment effect",
        title = "True treatment effect over study duration"
      ) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
      )
    
    p
    
  })
  
  output$outputPlot = renderPlot({
    
    data = simulated_data()
    
    p = data[["data"]] %>% 
      ggplot(aes(x = time, y = Yt)) + 
      geom_point() + 
      theme_minimal() +
      labs(
        x = "Observation number",
        y = "Observed outcome",
        title = "Observed outcome over study duration"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5)
      )
        
    p
    
  })
  
  
  
}
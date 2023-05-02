library(shiny)
library(tidyverse)

# Define UI for Shiny app
ui <- navbarPage("N-of-1 Data Generator",
                 tabPanel("Dataset Parameterization",
                          sidebarLayout(
                            sidebarPanel(
                              wellPanel(
                                h2("Trial Parameters"),
                                numericInput("n_trts", "Number of treatments", value = 2, min = 1),
                                numericInput("n_blocks", "Number of blocks", value = 3, min = 1),
                                numericInput("period_length", "Treatment period length", value = 10, min = 1),
                                numericInput("washout_length", "Washout period length", value = 5, min = 0),
                                numericInput("sampling_timestep", "Sampling time interval", value = 1, min = 0.1, step = 0.1)
                              ),
                              wellPanel(
                                h2("Treatment Parameters"),
                                h3("Treatment effects"),
                                uiOutput("Ejs"),
                                uiOutput("Gammas"),
                                uiOutput("Taus"),
                                numericInput("eta", "Time constant for the decay of the total treatment effect", value = 1, min = 0),
                                numericInput("sd_b", "Baseline noise", value = 1, min = 0),
                                numericInput("sd_p", "Process noise", value = 1, min = 0),
                                numericInput("sd_o", "Observation noise", value = 1, min = 0),
                                numericInput("phi", "Autocorrelation parameter (AR1)", value = 0.5)
                              ),
                              wellPanel(
                                h2("Simulation Parameters"),
                                numericInput("n_datasets", "Number of datasets to generate", value = 1)
                              ),
                              actionButton("simulate", "Simulate Datasets", width = "100%"),
                              actionButton("simulate_many", "Simulate Multiple Datasets", width = "100%")
                            ),
                            mainPanel(
                              h1("Dataset Plots"),
                              plotOutput("treatmentPlot"),
                              plotOutput("outputPlot"),
                              dataTableOutput("test")
                            )
                          )
                 )
)

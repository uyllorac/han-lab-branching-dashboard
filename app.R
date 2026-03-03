library(shiny)
library(ggplot2)
library(dplyr)
library(DT)

# Sample data based on your protocol (replace with your real CSV later)
sample_data <- data.frame(
  Branch_ID = paste0("B", 1:150),
  X = round(runif(150, 0, 240), 1),           # soma-to-tip position (pixels)
  Height = round(rnorm(150, 120, 40), 1),     # lifetime proxy (pixels/frames)
  Genotype = sample(c("WT", "PTP69D_KO"), 150, replace = TRUE, prob = c(0.5, 0.5)),
  Frequency = round(runif(150, 0.01, 0.15), 3)
)

ui <- fluidPage(
  titlePanel("Han Lab Dendritic Branching Dashboard - Interactive Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload your processed CSV (from protocol Step 3)", accept = ".csv"),
      selectInput("genotype", "Filter by Genotype:", 
                  choices = c("All", "WT", "PTP69D_KO"), selected = "All"),
      sliderInput("x_range", "Soma-to-Tip Position (X):", 
                  min = 0, max = 240, value = c(0, 240)),
      sliderInput("z_threshold", "Outlier Threshold (Z-score):", 
                  min = 1, max = 4, value = 2.5, step = 0.5),
      hr(),
      h5("Your Protocol Data"),
      p("Branch ID | X position | Height (lifetime) | Genotype | Frequency")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Interactive Plots",
                 plotOutput("violinPlot", height = "400px"),
                 plotOutput("scatterPlot", height = "400px")),
        tabPanel("Outliers & Summary",
                 DTOutput("outlierTable"),
                 verbatimTextOutput("summaryStats"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  data <- reactive({
    req(input$file)
    df <- read.csv(input$file$datapath)
    df$Genotype <- as.factor(df$Genotype)
    df
  }) %>% bindEvent(input$file)
  
  # Use sample data if no file uploaded yet
  current_data <- reactive({
    if (is.null(input$file)) sample_data else data()
  })
  
  filtered <- reactive({
    df <- current_data()
    if (input$genotype != "All") df <- df %>% filter(Genotype == input$genotype)
    df %>% filter(X >= input$x_range[1] & X <= input$x_range[2])
  })
  
  # Add Z-score for outlier detection (exactly like automated checks in J&J tools)
  with_outliers <- reactive({
    df <- filtered()
    df$z_score <- (df$Height - mean(df$Height, na.rm = TRUE)) / sd(df$Height, na.rm = TRUE)
    df$outlier <- abs(df$z_score) > input$z_threshold
    df
  })
  
  output$violinPlot <- renderPlot({
    ggplot(with_outliers(), aes(x = Genotype, y = Height, fill = Genotype)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_jitter(aes(color = outlier), width = 0.2, size = 2) +
      scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
      labs(title = "Branch Lifetime Distribution (Your Protocol Data)",
           y = "Height / Lifetime Proxy (pixels)") +
      theme_minimal()
  })
  
  output$scatterPlot <- renderPlot({
    ggplot(with_outliers(), aes(x = X, y = Height, color = Genotype)) +
      geom_point(size = 3, alpha = 0.8) +
      geom_smooth(method = "loess", se = FALSE) +
      labs(title = "Spatial Branching Frequency & Lifetime (Soma-to-Tip)",
           x = "Position from Soma (X pixels)") +
      theme_minimal()
  })
  
  output$outlierTable <- renderDT({
    with_outliers() %>% filter(outlier == TRUE) %>%
      select(Branch_ID, Genotype, X, Height, z_score) %>%
      datatable(options = list(pageLength = 10))
  })
  
  output$summaryStats <- renderPrint({
    df <- with_outliers()
    cat("Total branches:", nrow(df), "\n")
    cat("Flagged outliers:", sum(df$outlier), "\n")
    cat("Mean lifetime (Height):", round(mean(df$Height), 1), "\n")
    cat("WT mean:", round(mean(df$Height[df$Genotype=="WT"]), 1), "\n")
    cat("PTP69D_KO mean:", round(mean(df$Height[df$Genotype=="PTP69D_KO"]), 1))
  })
}

shinyApp(ui = ui, server = server)
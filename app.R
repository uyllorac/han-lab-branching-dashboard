library(shiny)
library(ggplot2)
library(dplyr)
library(DT)

sample_data <- data.frame(
  Branch_ID = paste0("B", 1:150),
  X = round(runif(150, 0, 240), 1),
  Height = round(rnorm(150, 120, 40), 1),
  Genotype = sample(c("WT", "PTP69D_KO"), 150, replace = TRUE, prob = c(0.5, 0.5)),
  Frequency = round(runif(150, 0.01, 0.15), 3)
)

sample_data_2nd <- data.frame(
  Branch_ID = paste0("B", 1:100),
  Genotype = sample(c("WT", "PTP69D_KO"), 100, replace = TRUE, prob = c(0.5, 0.5)),
  First_Emergence_Position = round(runif(100, 0, 240), 1),
  Lifetime = round(rnorm(100, 100, 35), 1),
  First_Emergence_Time = round(runif(100, 0, 500), 1)
)

ui <- fluidPage(
  titlePanel("Han Lab Dendritic Branching Dashboard - Interactive Analysis"),
  sidebarLayout(
    sidebarPanel(
      # --- Spatial file ---
      h4("Spatial Data"),
      fileInput("file", "Upload Spatial CSV (Step 3)", accept = ".csv"),
      selectInput("genotype", "Filter by Genotype:",
                  choices = c("All", "WT", "PTP69D_KO"), selected = "All"),
      sliderInput("x_range", "Soma-to-Tip Position (X):",
                  min = 0, max = 1200, value = c(0, 1200)),
      sliderInput("z_threshold", "Outlier Threshold (Z-score):",
                  min = 1, max = 4, value = 2.5, step = 0.5),
      hr(),
      # --- 2nd order file ---
      h4("2nd Order Branch Data"),
      fileInput("file_2nd", "Upload 2nd Order CSV", accept = ".csv"),
      selectInput("genotype_2nd", "Filter by Genotype (2nd Order):",
                  choices = c("All", "WT", "PTP69D_KO"), selected = "All"),
      sliderInput("z_threshold_2nd", "Outlier Threshold Z-score (2nd Order):",
                  min = 1, max = 4, value = 2.5, step = 0.5),
      hr(),
      h5("Your Protocol Data"),
      p("Branch ID | X position | Lifetime | Genotype")
    ),
    mainPanel(
      tabsetPanel(
        
        tabPanel("2nd Order Quantifications",
                 h4("Branch Lifetime — 2nd Order Branches"),
                 plotOutput("lifetimePlot_2nd", height = "400px"),
                 hr(),
                 h4("Branching Quantity per Genotype"),
                 plotOutput("quantityPlot_2nd", height = "400px")),
        
        tabPanel("Spatial Distribution",
                 sliderInput("norm_range", "Normalized X Range (%):",
                             min = 0, max = 100, value = c(0, 100), step = 5),
                 plotOutput("spatialJitter", height = "500px"),
                 plotOutput("spatialViolin", height = "500px")),
        
        tabPanel("Outliers & Summary",
                 DTOutput("outlierTable"),
                 verbatimTextOutput("summaryStats"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  # --- Spatial data pipeline ---
  data <- reactive({
    req(input$file)
    df <- read.csv(input$file$datapath)
    df$Genotype <- as.factor(df$Genotype)
    df
  }) %>% bindEvent(input$file)
  
  current_data <- reactive({
    if (is.null(input$file)) sample_data else data()
  })
  
  filtered <- reactive({
    df <- current_data()
    if (input$genotype != "All") df <- df %>% filter(Genotype == input$genotype)
    df %>% filter(X >= input$x_range[1] & X <= input$x_range[2])
  })
  
  with_outliers <- reactive({
    df <- filtered()
    df$z_score <- (df$Height - mean(df$Height, na.rm = TRUE)) / sd(df$Height, na.rm = TRUE)
    df$outlier <- abs(df$z_score) > input$z_threshold
    df
  })
  
  normalized <- reactive({
    df <- with_outliers()
    df <- df %>%
      group_by(Genotype) %>%
      mutate(X_norm = (X - min(X)) / (max(X) - min(X)) * 100) %>%
      ungroup()
    df %>% filter(X_norm >= input$norm_range[1] & X_norm <= input$norm_range[2])
  })
  
  # --- 2nd order data pipeline ---
  data_2nd <- reactive({
    req(input$file_2nd)
    df <- read.csv(input$file_2nd$datapath)
    df$Genotype <- as.factor(df$Genotype)
    # Normalize column names to remove spaces
    colnames(df) <- make.names(colnames(df))
    df
  }) %>% bindEvent(input$file_2nd)
  
  current_data_2nd <- reactive({
    if (is.null(input$file_2nd)) sample_data_2nd else data_2nd()
  })
  
  filtered_2nd <- reactive({
    df <- current_data_2nd()
    if (input$genotype_2nd != "All") df <- df %>% filter(Genotype == input$genotype_2nd)
    df
  })
  
  with_outliers_2nd <- reactive({
    df <- filtered_2nd()
    df$z_score <- (df$Lifetime - mean(df$Lifetime, na.rm = TRUE)) / sd(df$Lifetime, na.rm = TRUE)
    df$outlier <- abs(df$z_score) > input$z_threshold_2nd
    df
  })
  
  # --- 2nd order plots ---
  output$lifetimePlot_2nd <- renderPlot({
    ggplot(with_outliers_2nd(), aes(x = Genotype, y = Lifetime, fill = Genotype)) +
      geom_violin(trim = FALSE, alpha = 0.75) +
      geom_jitter(aes(color = outlier), width = 0.2, size = 2) +
      scale_fill_manual(values = c("WT" = "#F4A6A0", "PTP69D_KO" = "#7ECECA")) +
      scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
      labs(title = "2nd Order Branch Lifetime Distribution",
           y = "Lifetime (pixels)") +
      theme_minimal()
  })
  
  output$quantityPlot_2nd <- renderPlot({
    count_df <- with_outliers_2nd() %>%
      group_by(Genotype) %>%
      summarise(Count = n(), .groups = "drop")
    
    ggplot(count_df, aes(x = Genotype, y = Count, fill = Genotype)) +
      geom_bar(stat = "identity", alpha = 0.85, width = 0.5) +
      geom_text(aes(label = Count), vjust = -0.5, size = 5) +
      scale_fill_manual(values = c("WT" = "#F4A6A0", "PTP69D_KO" = "#7ECECA")) +
      labs(title = "2nd Order Branching Quantity per Genotype",
           y = "Number of Branches") +
      theme_minimal()
  })
  
  # --- Spatial plots ---
  output$spatialJitter <- renderPlot({
    ggplot(normalized(), aes(x = Genotype, y = X_norm, color = Genotype)) +
      geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
      scale_y_continuous(
        labels = function(x) paste0(x, "%"),
        limits = c(input$norm_range[1], input$norm_range[2])
      ) +
      scale_color_manual(values = c("WT" = "#F4A6A0", "PTP69D_KO" = "#7ECECA")) +
      labs(title = "Spatial Distribution (Normalized 0-100% within Group)",
           subtitle = "WT vs PTP69D_KO",
           x = "Group", y = "Normalized X (0-100%)", color = "Group") +
      theme_minimal()
  })
  
  output$spatialViolin <- renderPlot({
    ggplot(normalized(), aes(x = Genotype, y = X_norm, fill = Genotype)) +
      geom_violin(trim = FALSE, alpha = 0.85) +
      geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
      stat_summary(fun = mean, geom = "point", color = "red", size = 3) +
      scale_y_continuous(
        labels = function(x) paste0(x, "%"),
        limits = c(input$norm_range[1], input$norm_range[2])
      ) +
      scale_fill_manual(values = c("WT" = "#F4A6A0", "PTP69D_KO" = "#7ECECA")) +
      labs(title = "Spatial Distribution (Normalized 0-100% within Group)",
           subtitle = "WT vs PTP69D_KO",
           x = "Group", y = "Normalized X (0-100%)", fill = "Group") +
      theme_minimal()
  })
  
  # --- Outliers & Summary (spatial data) ---
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

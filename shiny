# Shiny App 
library(shiny)
library(ggplot2)
library(reshape2)
library(dplyr)

# FUNCTIONS (loop, vectorized, modular versions)

# ------- LOOP VERSION -------
# Computes the similarity matrix using fully iterative loops
similarity_loop <- function(dat) {
  n <- nrow(dat)
  S <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      xi <- dat[i, ]; xj <- dat[j, ]
      xi[is.na(xi)] <- mean(xi, na.rm = TRUE)  
      xj[is.na(xj)] <- mean(xj, na.rm = TRUE)
      Dist    <- sqrt(sum((xi - xj)^2))
      Corr    <- cor(xi, xj)
      MeanDif <- mean(abs(xi - xj))
      Thresh  <- (mean(xi) + mean(xj)) / 2
      Count   <- sum((xi > Thresh) & (xj > Thresh))
      S[i, j] <- (Dist + abs(Corr) + MeanDif) / (Count + 1)
    }
  }
  S
}

# ------- VECTORIZED VERSION -------
# Faster implementation using vectorized operations and array expansion
similarity_vectorized <- function(dat) {
  dat_imp <- apply(dat, 2, function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x })
  n <- nrow(dat_imp); p <- ncol(dat_imp)
  means <- rowMeans(dat_imp)
  
  Dist    <- as.matrix(dist(dat_imp))           
  Corr    <- abs(cor(t(dat_imp)))             
  MeanDif <- as.matrix(dist(dat_imp, method = "manhattan")) / p
  
  Thresh_mat <- outer(means, means, FUN = function(a, b) (a + b)/2)
  
  A <- array(rep(dat_imp, times = n), dim = c(n, p, n))
  B <- aperm(A, c(3,2,1))
  Thresh3D <- array(rep(Thresh_mat, each = p), dim = c(n, p, n))
  
  Above <- (A > Thresh3D) & (B > Thresh3D)
  Count <- apply(Above, c(1,3), sum, na.rm = TRUE)
  
  (Dist + Corr + MeanDif) / (Count + 1)
}

# ------- MODULAR VERSION -------
# Helper functions for a more readable modular structure
eucl <- function(a, b) sqrt(sum((a - b)^2))
corr_mod <- function(a, b) abs(cor(a, b))
mean_abs <- function(a, b) mean(abs(a - b))
count_thr <- function(a, b) { t <- (mean(a) + mean(b)) / 2; sum((a > t) & (b > t)) }

# Combines the modular components into the final similarity matrix
similarity_modular <- function(dat) {
  n <- nrow(dat)
  dat_imp <- apply(dat, 2, function(x){ x[is.na(x)] <- mean(x, na.rm=TRUE); x })
  S <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      S[i,j] <- (eucl(dat_imp[i,], dat_imp[j,]) +
                   corr_mod(dat_imp[i,], dat_imp[j,]) +
                   mean_abs(dat_imp[i,], dat_imp[j,])) /
        (count_thr(dat_imp[i,], dat_imp[j,]) + 1)
    }
  }
  S
}


# --------------------------
# UI
# --------------------------
ui <- fluidPage(
  titlePanel("Similarity Index Calculator – Shiny App"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload CSV Data", accept = ".csv"),
      
      selectInput("version", "Choose Method:",
                  choices = c("Loop Version", "Vectorized Version", "Modular Version")),
      
      h4("Skewed markers"),
      checkboxGroupInput("skewed_vars", NULL,
                         choices = c("HER2","Tumor_size","p53","CEA","CA125","CA199",
                                     "AFP","hCG","LDH","ALP","Ferritin","CRP","IL6")),
      
      h4("Proportional markers"),
      checkboxGroupInput("prop_vars", NULL,
                         choices = c("Ki67","ctDNA")),
      
      h4("Tumor density"),
      checkboxGroupInput("density_vars", NULL,
                         choices = c("Tumor_density")),
      
      h4("Tumor growth"),
      checkboxGroupInput("growth_vars", NULL,
                         choices = c("Tumor_growth")),
      
      h4("ESR"),
      checkboxGroupInput("esr_vars", NULL,
                         choices = c("ESR")),
      
      # Quick selection buttons
      actionButton("select_all", "Select All Variables"),
      actionButton("deselect_all", "Deselect All Variables"),
      
      sliderInput("n_rows", "Number of rows to preview:",
                  min = 5, max = 90, value = 10),
      
      checkboxInput("standardize", "Standardize data before analysis", value = TRUE),
      
      actionButton("run", "Compute Similarity Matrix"),
      
      br(), br(),
      helpText("Select variables from each group to compute similarity indices.")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Heatmap", plotOutput("heatmap")),
        tabPanel("Execution Time", tableOutput("time_table")),
        tabPanel("Data Preview",
                 tableOutput("preview"),
                 h4("Summary of Selected Variables"),
                 tableOutput("summary"))
      )
      
    )
  )
)


# --------------------------
# SERVER
# --------------------------
server <- function(input, output, session) {
  
  # Load uploaded data
  data_input <- reactive({
    req(input$file)
    read.csv(input$file$datapath, check.names = TRUE)
  })
  
  # Data preview
  output$preview <- renderTable({
    head(data_input(), input$n_rows)
  })
  
  # Summary table of selected variables
  output$summary <- renderTable({
    df <- data_input()
    df <- df[, -1, drop = FALSE]   # remove ID column
    
    selected_vars <- c(input$skewed_vars,
                       input$prop_vars,
                       input$density_vars,
                       input$growth_vars,
                       input$esr_vars)
    
    selected_vars <- intersect(names(df), selected_vars)
    if (length(selected_vars) == 0) return(NULL)
    
    dat <- head(df[, selected_vars, drop = FALSE], input$n_rows)
    
    # Basic descriptive statistics
    summary_df <- data.frame(
      Variable = names(dat),
      Mean = sapply(dat, function(x) round(mean(x, na.rm = TRUE), 2)),
      SD   = sapply(dat, function(x) round(sd(x, na.rm = TRUE), 2)),
      Min  = sapply(dat, function(x) round(min(x, na.rm = TRUE), 2)),
      Max  = sapply(dat, function(x) round(max(x, na.rm = TRUE), 2))
    )
    summary_df
  })
  
  # "Select all" buttons
  observeEvent(input$select_all, {
    updateCheckboxGroupInput(session, "skewed_vars",
                             selected = c("HER2","Tumor_size","p53","CEA","CA125","CA199",
                                          "AFP","hCG","LDH","ALP","Ferritin","CRP","IL6"))
    updateCheckboxGroupInput(session, "prop_vars", selected = c("Ki67","ctDNA"))
    updateCheckboxGroupInput(session, "density_vars", selected = c("Tumor_density"))
    updateCheckboxGroupInput(session, "growth_vars", selected = c("Tumor_growth"))
    updateCheckboxGroupInput(session, "esr_vars", selected = c("ESR"))
  })
  
  observeEvent(input$deselect_all, {
    updateCheckboxGroupInput(session, "skewed_vars", selected = character(0))
    updateCheckboxGroupInput(session, "prop_vars", selected = character(0))
    updateCheckboxGroupInput(session, "density_vars", selected = character(0))
    updateCheckboxGroupInput(session, "growth_vars", selected = character(0))
    updateCheckboxGroupInput(session, "esr_vars", selected = character(0))
  })
  
  observeEvent(input$run, {
    df <- data_input()
    df <- df[, -1, drop = FALSE]  # remove ID column
    
    selected_vars <- c(input$skewed_vars,
                       input$prop_vars,
                       input$density_vars,
                       input$growth_vars,
                       input$esr_vars)
    
    selected_vars <- intersect(names(df), selected_vars)
    
    # Handle empty selection
    if (length(selected_vars) == 0) {
      output$time_table <- renderTable(data.frame(Message="No variables selected"))
      output$heatmap <- renderPlot({ plot.new(); text(0.5, 0.5, "No variables selected") })
      return(NULL)
    }
    
    dat <- head(df[, selected_vars, drop = FALSE], input$n_rows)
    dat <- as.matrix(dat)
    storage.mode(dat) <- "double"
    
    if (input$standardize) dat <- scale(dat)
    
    # Benchmark execution time of each method
    times <- data.frame(Method = character(), Time = numeric(), stringsAsFactors = FALSE)
    
    start <- Sys.time(); similarity_loop(dat); end <- Sys.time()
    times <- rbind(times, data.frame(Method = "Loop Version", 
                                     Time = as.numeric(difftime(end, start, units="secs"))))
    
    start <- Sys.time(); similarity_vectorized(dat); end <- Sys.time()
    times <- rbind(times, data.frame(Method = "Vectorized Version", 
                                     Time = as.numeric(difftime(end, start, units="secs"))))
    
    start <- Sys.time(); similarity_modular(dat); end <- Sys.time()
    times <- rbind(times, data.frame(Method = "Modular Version", 
                                     Time = as.numeric(difftime(end, start, units="secs"))))
    
    # Output timing table
    output$time_table <- renderTable({
      times
    })
    
    # Compute similarity only for the chosen method 
    S <- switch(input$version,
                "Loop Version" = similarity_loop(dat),
                "Vectorized Version" = similarity_vectorized(dat),
                "Modular Version" = similarity_modular(dat))
    
    # Heatmap
    output$heatmap <- renderPlot({
      S_melt <- melt(S)
      ggplot(S_melt, aes(Var1, Var2, fill = value)) +
        geom_tile() +
        scale_fill_viridis_c() +
        theme_minimal() +
        labs(title = paste("Heatmap - Variables:", paste(selected_vars, collapse=", ")))
    })
  })
}

shinyApp(ui, server)

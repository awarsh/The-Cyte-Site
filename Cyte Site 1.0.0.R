#This is an app to showcase data

library(shiny) 
library(shinythemes)
library(readxl) 
library(ggplot2)
library(plotly)
library(ggpattern)
library(readr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(janitor)
library(patchwork)
library(shinydashboard)
library(DT)
library(scales)


# Now that we have our tools ready, lets bring something in to work on.
data <- read_csv("~/Downloads/hannahs_processed_data_no_NA_NaN.csv")

stats <- read_csv("~/Downloads/Hannahs_stats.csv") %>%
  janitor::clean_names()

# Create a shinydashboard UI
ui <- dashboardPage(
  dashboardHeader(title = "Cyte Site"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Project Description", tabName = "project_desc", icon = icon("info-circle")),
      menuItem("Data Analysis", tabName = "data_analysis", icon = icon("chart-bar"))
    )
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "project_desc",
              h3("Project Description"),
              p("This interactive application is designed to explore comprehensive proteomic data across four different cell types: Pericytes, BMVEC, Hepatocytes, and Astrocytes. Dive into the world of proteins, genes, and their varying expressions in different cell environments.")
      ),
      
      tabItem(tabName = "data_analysis",
              fluidRow(
                column( width = 12,
                        box(
                          title = "Input Selection",
                          status = "primary",
                          radioButtons("selection", "Select an option:",
                                       choices = list("Gene" = "gene", "Protein" = "protein")),
                          uiOutput("dynamic_ui"),
                          checkboxGroupInput("cell_choice", "Select cell types of interest:",
                                             choices = c("BMVECs", "Pericytes", "Hepatocytes", "Astrocytes"),
                                             selected = c("BMVECs", "Pericytes", "Hepatocytes", "Astrocytes")),
                        ),
                        box(
                          title = "Protein Information",
                          status = "primary",
                          solidHeader = TRUE,
                          collapsible = TRUE,
                          htmlOutput("proteinInfo")  
                        ),
                        box(
                          title = "MS1 Statistics",
                          status = "primary",
                          solidHeader = TRUE,
                          collapsible = TRUE,
                          uiOutput("statistics")  
                        )
                )
              ),
              
              fluidRow(
                box(
                  title = "Data Analysis",
                  width = 12,  
                  tabsetPanel(
                    tabPanel("Plots", 
                             radioButtons("log2_choice_bar", "Choose Scale:",
                                          choices = c("Original" = "original", "Log2" = "log2")),
                             fluidRow(
                               column(4, plotlyOutput("ms1_plot")),
                               column(4, plotlyOutput("concentration_plot")),
                               column(4, plotlyOutput("copy_number_plot"))
                             ),
                             checkboxInput("show_points", "Show individual points")
                             
                    )
                  )
                  
                )
              )
      )
    )
  )
)


server <- function(input, output) {
  
  # Render UI for dynamic selection of gene or protein
  output$dynamic_ui <- renderUI({
    if (input$selection == "gene") {
      selectInput("gene_choice", "Gene", choices = unique(data$pg_genes))
    } else {
      selectInput("protein_choice", "Protein", choices = unique(data$pg_uni_prot_ids))
    }
  })
  
  # Observe and print the selection
  observeEvent(input$gene_choice, {
    print(input$gene_choice)
  })
  
  observeEvent(input$protein_choice, {
    print(input$protein_choice)
  })
  
  #filter data based on the ui selections
  filtered_data <- reactive({
    
    selected_cell_type <- input$cell_choice
    selected_gene <- input$gene_choice
    selected_protein <- input$protein_choice
    
    filtered_data <- data %>%
      filter(cell_type %in% selected_cell_type)
    
    if(input$selection == "gene") {
      filtered_data <- filtered_data %>% filter(pg_genes %in% selected_gene)
    } else {
      filtered_data <- filtered_data %>% filter(pg_uni_prot_ids %in% selected_protein)
    }
    
    filtered_data
  })
  
  #create the summarized data to be used in the graphs
  summary_data <- reactive({
    selected_column <- if (input$selection == "gene") "pg_genes" else "pg_protein_groups"
    
    # Summarize the data
    filtered_data() %>%
      mutate(selected_col_value = .data[[selected_column]]) %>%
      group_by(cell_type, data_type, selected_col_value) %>%
      mutate(mean_value = mean(value), sd_value = sd(value), n = n()) %>%
      mutate(se = sd_value / sqrt(n)) %>%
      ungroup()
  })
  
  output$ms1_plot <- renderPlotly({
    legend_title <- if (input$selection == "gene") "Gene" else "Protein"
    
    # If log2 transformation is selected, transform the individual data points first.
    if (input$log2_choice_bar == "log2") {
      plot_data <- summary_data() %>%
        mutate(value = log2(value + 1)) # Transform individual data points
    } else {
      plot_data <- summary_data()
    }
    
    # After transformation, recalculate mean_value, se, and error bars on the transformed scale.
    plot_data <- plot_data %>%
      filter(data_type == "MS1 Quantity") %>%
      group_by(cell_type, data_type, selected_col_value) %>%
      mutate(
        mean_value = mean(value), # Recalculate mean_value on the transformed scale
        se = sd(value) / sqrt(n())  # Recalculate se on the transformed scale
      ) %>%
      ungroup() %>%
      mutate(
        lower = pmax(mean_value - se, 0), # Ensure lower bound is non-negative
        upper = mean_value + se
      )
    
    p <- ggplot(plot_data, aes(x = cell_type, fill = cell_type, y = mean_value)) + 
      geom_bar(stat = "identity", position = "dodge") +
      geom_errorbar(aes(ymin = lower, ymax = upper),  # Updated aesthetics
                    position = position_dodge(0.9),
                    width = 0.25) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10))) +
      scale_fill_brewer(name = "Cell Types", palette = "Pastel1") +
      scale_y_continuous(labels = label_comma()) +
      labs(y = "MS1 Quantity", x = "") +
      facet_wrap(~selected_col_value)
    
    if (input$show_points) {
      p <- p + geom_point(data = plot_data, aes(x = cell_type, y = value), 
                          position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9))
    }
    
    plotly::ggplotly(p)
  })
  
  output$concentration_plot <- renderPlotly({
    legend_title <- if (input$selection == "gene") "Gene" else "Protein"
    
    # If log2 transformation is selected, transform the individual data points first.
    if (input$log2_choice_bar == "log2") {
      plot_data <- summary_data() %>%
        mutate(value = log2(value + 1)) # Transform individual data points
    } else {
      plot_data <- summary_data()
    }
    
    # After transformation, recalculate mean_value, se, and error bars on the transformed scale.
    plot_data <- plot_data %>%
      filter(data_type == "Concentration") %>%
      group_by(cell_type, data_type, selected_col_value) %>%
      mutate(
        mean_value = mean(value), # Recalculate mean_value on the transformed scale
        se = sd(value) / sqrt(n())  # Recalculate se on the transformed scale
      ) %>%
      ungroup() %>%
      mutate(
        lower = pmax(mean_value - se, 0), # Ensure lower bound is non-negative
        upper = mean_value + se
      )
    
    p <- ggplot(plot_data, aes(x = cell_type, fill = cell_type, y = mean_value)) + 
      geom_bar(stat = "identity", position = "dodge") +
      geom_errorbar(aes(ymin = lower, ymax = upper),  # Updated aesthetics
                    position = position_dodge(0.9),
                    width = 0.25) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10))) +
      scale_fill_brewer(name = "Cell Types", palette = "Pastel1") +
      scale_y_continuous(labels = label_comma()) +
      labs(y = "Concentration (nM)", x = "") +
      facet_wrap(~selected_col_value)
    
    if (input$show_points) {
      p <- p + geom_point(data = plot_data, aes(x = cell_type, y = value), 
                          position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9))
    }
    
    plotly::ggplotly(p)
  })
  
  output$copy_number_plot <- renderPlotly({
    legend_title <- if (input$selection == "gene") "Gene" else "Protein"
    
    # If log2 transformation is selected, transform the individual data points first.
    if (input$log2_choice_bar == "log2") {
      plot_data <- summary_data() %>%
        mutate(value = log2(value + 1)) # Transform individual data points
    } else {
      plot_data <- summary_data()
    }
    
    # After transformation, recalculate mean_value, se, and error bars on the transformed scale.
    plot_data <- plot_data %>%
      filter(data_type == "Copy Number") %>%
      group_by(cell_type, data_type, selected_col_value) %>%
      mutate(
        mean_value = mean(value), # Recalculate mean_value on the transformed scale
        se = sd(value) / sqrt(n())  # Recalculate se on the transformed scale
      ) %>%
      ungroup() %>%
      mutate(
        lower = pmax(mean_value - se, 0), # Ensure lower bound is non-negative
        upper = mean_value + se
      )
    
    p <- ggplot(plot_data, aes(x = cell_type, fill = cell_type, y = mean_value)) + 
      geom_bar(stat = "identity", position = "dodge") +
      geom_errorbar(aes(ymin = lower, ymax = upper),  # Updated aesthetics
                    position = position_dodge(0.9),
                    width = 0.25) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10))) +
      scale_fill_brewer(name = "Cell Types", palette = "Pastel1") +
      scale_y_continuous(labels = label_comma()) +
      labs(y = "Copy Number (copies/cell)", x = "") +
      facet_wrap(~selected_col_value)
    
    if (input$show_points) {
      p <- p + geom_point(data = plot_data, aes(x = cell_type, y = value), 
                          position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9))
    }
    
    plotly::ggplotly(p)
  })
  
  # Reactive expression for the table data
  table_data <- reactive({
    # Start with the full dataset
    temp_data <- data  # make sure 'data' is available in the server environment
    
    # Filter by selected gene or protein
    if (input$selection == "gene") {
      selected_gene <- input$gene_choice
      temp_data <- temp_data %>% filter(pg_genes == selected_gene)
    } else {
      selected_protein <- input$protein_choice
      temp_data <- temp_data %>% filter(pg_protein_groups == selected_protein)
    }
    
    # Select only the relevant columns
    temp_data <- temp_data %>% 
      select(pg_genes, pg_uni_prot_ids, pg_protein_names, pg_molecular_weight, pg_protein_descriptions, histones) %>% 
      slice(1)
    
    temp_data  # Return the filtered data
  })
  
  #create the protein infor table
  output$proteinInfo <- renderUI({
    data <- table_data()
    
    tags$dl(
      tags$dd(strong("Gene:"), data$pg_genes),
      tags$dd(strong("Protein ID:"), data$pg_uni_prot_ids),
      tags$dd(strong("Protein Name:"), data$pg_protein_names),
      tags$dd(strong("Molecular Weight:"), data$pg_molecular_weight, "g/mol"),
      tags$dd(strong("Description:"), data$pg_protein_descriptions)
    )
  })
  
  #filter the statistics based on the user inputs 
  filtered_stats <- reactive({
    if (input$selection == "gene") {
      stats %>% 
        filter(genes == input$gene_choice) %>%
        select(comparison_group1_group2, pvalue, qvalue)
    } else {
      stats %>% 
        filter(protein_groups == input$protein_choice) %>%
        select(comparison_group1_group2, pvalue, qvalue)
    }
  })
  
  #output the statistics
  output$statistics <- renderUI({
    # Call the reactive expression to get the filtered data
    stats_data <- filtered_stats()
    
    if (nrow(stats_data) > 0) {
      # Create a div for each element and use inline-block for horizontal layout
      stats_divs <- lapply(1:nrow(stats_data), function(i) {
        tags$div(
          style = "display: inline-block; vertical-align: top; margin-right: 20px;",
          tags$div(strong("Comparison: "), stats_data[[i, "comparison_group1_group2"]]),
          tags$div(strong("P-value: "), format(stats_data[[i, "pvalue"]], scientific = TRUE)),
          tags$div(strong("Q-value: "), format(stats_data[[i, "qvalue"]], scientific = TRUE))
        )
      })
      
      # Wrap the divs in a container with horizontal scrolling
      tags$div(
        style = "width: 100%; overflow-x: auto; white-space: nowrap;",
        do.call(tagList, stats_divs)
      )
    } else {
      tags$p("No statistical data available for this selection.")
    }
  })
}
# Run the app
shinyApp(ui = ui, server = server)

library(shinydashboard)
library(rhandsontable)
library(shinyjs)
library(shiny)
library(CytoNorm)
library(flowCore)
library(DT)
library(shinycssloaders)

# =============================================
# Dashboard UI for cyCombine
# =============================================
ui <- dashboardPage(
  skin = "blue",
  
  # ---------------- Header ----------------
  dashboardHeader(title = "cyCombine"),
  
  # ---------------- Sidebar ----------------
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data", tabName = "data", icon = icon("folder-open")),
      menuItem("Run cyCombine", tabName = "run", icon = icon("play"))
    )
  ),
  
  # ---------------- Body ----------------
  dashboardBody(
    tabItems(
      
      # ================= DATA TAB =================
      tabItem(
        tabName = "data",
        fluidPage(
          
          # ---------------- Upload & Experiment Info ----------------
          box(
            width = 10,
            
            tags$h3(
              "Upload Data",
              style = "font-weight: bold; background-color: #fff9c4; padding: 4px; display: inline-block;"
            ),
            
            # Experiment name (default: timestamp)
            textInput(
              "experiment_name",
              "Experiment name",
              value = format(Sys.time(), "%Y%m%d_%H%M%S")
            ),
            
            # Number of batches
            numericInput(
              "n_batches",
              "Number of batches",
              value = 2,
              min = 1,
              width = "20%"
            ),
            
            # Dynamic UI for batch uploads
            uiOutput("batch_upload_ui"),
            
            tags$br(),
            tags$br(),
            
            # ---------------- Panel Harmonization ----------------
            tags$h3(
              "Panel Harmonization",
              style = "font-weight: bold; background-color: #fff9c4; padding: 4px; display: inline-block;"
            ),
            tags$br(),
            
            downloadButton("export_panel", "Export panel (Excel)"),
            fileInput("import_panel", "Import modified panel", accept = ".xlsx"),
            
            # Editable panel table
            rHandsontableOutput("panel_table", width = "100%", height = "400px"),
            tags$br(),
            tags$br(),
            
            tags$h5("Please click the 'Validate Design' button before proceeding"),
            uiOutput("covar_ui"),
            
            actionButton(
              "Validate_design",
              "Validate Design",
              style = "background-color: #FF5733; color: white; font-weight: bold; font-size: 16px; padding: 10px 20px; border-radius: 5px;"
            )
          ),
          
          # ---------------- Transformation Box ----------------
          box(
            width = 10,
            fluidRow(
              
              # Left column: transformation heading
              column(
                width = 10,
                tags$h3(
                  "Transformation",
                  style = "font-weight: bold; background-color: #fff9c4; padding: 4px; display: inline-block;"
                ),
                br()
              )
            ),
            
            fluidRow(
              
              # Right column: Transformation mode & cofactor
              column(
                width = 3,
                
                # Transformation mode: global or batch-specific
                radioButtons(
                  "transform_mode",
                  "Transformation mode",
                  choices = c(
                    "Same cofactor for all batches" = "global",
                    "Different cofactor per batch" = "batch"
                  ),
                  selected = "global"
                ),
                
                # Global cofactor input
                conditionalPanel(
                  condition = "input.transform_mode == 'global'",
                  numericInput("global_cofactor", "Cofactor", value = 500, min = 0.1)
                ),
                
                # Batch-specific cofactors (dynamic UI)
                conditionalPanel(
                  condition = "input.transform_mode == 'batch'",
                  uiOutput("batch_cofactor_ui")
                )
              ),
              
              # Middle column: Apply transformation button
              column(
                width = 6,
                actionButton(
                  "apply_transform",
                  "Apply transformation"
                )
              )
            )
          )
        )
      ),
      
      # ================= RUN TAB =================
      tabItem(
        tabName = "run",
        fluidPage(
          
          # ---------------- Run cyCombine ----------------
          box(
            width = 12,
            tags$h3(
              "Run CyCombine",
              style = "font-weight: bold; background-color: #fff9c4; padding: 4px; display: inline-block;"
            ),
            
            # UI for selecting reference batch
            uiOutput("goal_ref_batch_ui"),
            tags$br(),
            
            # Normalization method: scale or rank
            radioButtons(
              inputId = "norm_method",
              label = "Norm method. *Rank is recommended for heavy batch effects",
              choices = c("scale", "rank"),
              inline = TRUE
            ),
            
            # SOM iterations
            numericInput(
              "rlen",
              "Number of iterations for SOM training (*Higher values may improve results)",
              value = 10,
              min = 10
            ),
            
            # Run button
            actionButton("run_cyCombine", "Run cyCombine"),
            br(), br(),
            
            # ---------------- Export normalized files ----------------
            tags$h3(
              "Export normalized files",
              style = "font-weight: bold; background-color: #fff9c4; padding: 4px; display: inline-block;"
            ),
            tags$br(),
            downloadButton("download_norm", "Download normalized files")
          ),
          
          # ---------------- Visualization & Plots ----------------
          uiOutput("visualization_ui"),
          uiOutput("dens_plot_ui")
        )
      )
    )
  )
)
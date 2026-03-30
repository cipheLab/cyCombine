library(shiny)
library(cyCombine)
library(flowCore)
library(DT)
library(shinycssloaders)
library(openxlsx)
library(shinydashboard)
library(dplyr)
library(ggplot2)
library(flowDensity)

options(expressions = 5e5, shiny.maxRequestSize = 100 * 1024^3)

server <- function(input, output, session) {
  # ---------------------------------------------
  # Function to get batch name
  # Returns the user-provided name, or "Batch i" if empty
  # ---------------------------------------------
  get_batch_name <- function(i) {
    name <- rv$batch_names[[i]]
    if (is.null(name) || name == "") {
      return(paste("Batch", i))
    }
    return(name)
  }
  
  # ---------------------------------------------
  # Directory to save models
  # Creates directory if it doesn't exist
  # ---------------------------------------------
  model_dir <- "/media/data/html/Database/cytoNorm_Models"
  
  if (!dir.exists(model_dir)) {
    dir.create(model_dir, recursive = TRUE)
  }
  
  # ---------------------------------------------
  # Reactive values initialization
  # Stores uploaded data, design table, logs, transformations, batch names, etc.
  # ---------------------------------------------
  rv <- reactiveValues(
    fs_list = list(),          # List of flowSet objects per batch
    design = data.frame(),     # Design table with batch info
    log = "",                  # Log messages
    transfoTable = NULL,       # Table of marker transformations
    batch_names = list(),      # Names of batches
    isTransformed = FALSE,     # Flag for whether data is transformed
    transformed_markers = NULL # Markers that have been transformed
  )
  
  # ---------------------------------------------
  # Render UI for batch uploads dynamically
  # Depends on number of batches specified by user
  # ---------------------------------------------
  output$batch_upload_ui <- renderUI({
    req(input$n_batches)
    
    tagList(
      lapply(1:input$n_batches, function(i) {
        box(
          title = paste("Batch", i),
          width = 12,
          
          # Text input for batch name
          textInput(
            paste0("batch_name_", i),
            label = "Batch name",
            value = paste0("Batch", i)
          ),
          
          # File input for FCS files
          tags$div(
            style = "background-color: #f9f9f9; padding: 8px; border-radius: 5px;",
            fileInput(
              inputId = paste0("batch_", i),
              label = paste("Upload files for batch", i),
              multiple = TRUE,
              accept = ".fcs",
              buttonLabel = "Choose FCS files",
              placeholder = "No files selected",
              width = "100%"
            )
          ),
          
          # Placeholder for additional UI (e.g., reference selection)
          uiOutput(paste0("ref_ui_", i))
        )
      })
    )
  })
  
  # ---------------------------------------------
  # Update batch names reactive value whenever input changes
  # ---------------------------------------------
  observe({
    req(input$n_batches)
    
    for (i in 1:input$n_batches) {
      rv$batch_names[[i]] <- input[[paste0("batch_name_", i)]]
    }
  })
  
  # ---------------------------------------------
  # Observe file uploads for each batch
  # Reads FCS files, applies compensation, updates design table
  # ---------------------------------------------
  observe({
    req(input$n_batches)
    
    for (i in 1:input$n_batches) {
      local({
        batch_i <- i
        observeEvent(input[[paste0("batch_", batch_i)]], {
          files_input <- input[[paste0("batch_", batch_i)]]
          req(files_input)
          
          # Read each FCS file
          fs_list <- lapply(seq_len(nrow(files_input)), function(j) {
            ff <- read.FCS(
              files_input$datapath[j],
              transformation = FALSE,
              truncate_max_range = FALSE
            )
            
            # -------------------------
            # Automatic compensation
            # -------------------------
            spill <- FlowCIPHE::found.spill.CIPHE(ff)
            
            if (!is.null(spill) && length(spill) > 0 && !all(is.na(spill)) && spill[1] != "NULL") {
              spill_name <- spill[1]
              
              if (spill_name %in% names(ff@description)) {
                print(paste("Applying compensation with:", spill_name))
                ff <- flowCore::compensate(ff, ff@description[[spill_name]])
              } else {
                print("Spill name not found in description")
              }
              
            } else {
              print("No valid spill matrix found")
            }
            
            return(ff)
          })
          
          # Store processed FCS files in reactive values
          rv$fs_list[[batch_i]] <- fs_list
          
          # Update design table
          rv$design <- rv$design[rv$design$Batch != batch_i, ]
          new_rows <- data.frame(
            Batch = batch_i,
            BatchName = get_batch_name(batch_i),
            File = files_input$name,
            Type = "Sample",
            stringsAsFactors = FALSE
          )
          rv$design <- rbind(rv$design, new_rows)
          rownames(rv$design) <- NULL
          
          # Initialize transformation table based on first batch
          markers <- colnames(exprs(rv$fs_list[[1]][[1]]))
          names(markers) <- NULL
          rv$transfoTable <- data.frame(
            Fluo = markers,
            Arg = rep("none", length(markers))
          )
          
        }, ignoreNULL = TRUE)
      })
    }
  })
  
  # ---------------------------------------------
  # Build panel table for all batches
  # Shows fluorochrome and marker names for comparison
  # ---------------------------------------------
  observe({
    req(length(rv$fs_list) > 0)
    
    n_batches <- length(rv$fs_list)
    panel_list <- list()
    max_len <- 0
    
    for (i in seq_len(n_batches)) {
      ff <- rv$fs_list[[i]][[1]]
      param_data <- ff@parameters@data
      
      fluo <- param_data$name
      marker <- param_data$desc
      
      # Replace empty marker names with fluorochrome names
      marker[is.na(marker) | marker == ""] <- fluo
      
      # Detect if using mass cytometry (e.g., "123X_...")
      is_mass <- mean(grepl("^[0-9]+[A-Za-z]+_", marker)) > 0.1
      if (is_mass) {
        marker <- sub(".*_", "", marker)
      }
      
      df <- data.frame(
        Fluo = fluo,
        Marker = marker,
        stringsAsFactors = FALSE
      )
      
      panel_list[[i]] <- df
      max_len <- max(max_len, nrow(df))
    }
    
    # Merge batch panels into one table
    panel_table <- data.frame(matrix(nrow = max_len, ncol = 0))
    
    for (i in seq_len(n_batches)) {
      df <- panel_list[[i]]
      
      if (nrow(df) < max_len) {
        df[(nrow(df)+1):max_len, ] <- NA
      }
      
      colnames(df) <- c(
        paste0(get_batch_name(i), "_Fluo"),
        paste0(get_batch_name(i), "_Marker")
      )
      
      panel_table <- cbind(panel_table, df)
    }
    
    panel_table$Normalize <- FALSE
    rv$panel_table <- panel_table
  })

  # ---------------------------------------------
  # Render interactive panel table with rhandsontable
  # Allows editing directly in the Shiny app
  # ---------------------------------------------
  output$panel_table <- renderRHandsontable({
    req(rv$panel_table)
    
    rhandsontable(rv$panel_table, stretchH = "all") %>%
      hot_table(highlightCol = TRUE, highlightRow = TRUE)
  })
  
  # Update reactive panel table when user edits it
  observeEvent(input$panel_table, {
    req(input$panel_table)
    rv$panel_table <- hot_to_r(input$panel_table)
  })
  
  # ---------------------------------------------
  # Update sample/reference type based on user selection
  # Loops through batches and updates rv$design$Type
  # ---------------------------------------------
  observe({
    req(input$n_batches)
    for (i in 1:input$n_batches) {
      local({
        batch_i <- i
        ref_input_id <- paste0("ref_files_", batch_i)
        selected_ref <- input[[ref_input_id]]
        if (!is.null(selected_ref)) {
          idx <- rv$design$Batch == batch_i
          rv$design$Type[idx] <- ifelse(
            rv$design$File[idx] %in% selected_ref,
            "Reference",
            "Sample"
          )
        }
      })
    }
  })
  
  # ---------------------------------------------
  # Export transformation table to Excel
  # ---------------------------------------------
  output$exportTransformation <- downloadHandler(
    filename = function() {
      paste("transf_table.xlsx", sep = "")
    },
    content = function(file) {
      write.xlsx(rv$transfoTable, file)
    }
  )
  
  # Update transformation table from rhandsontable input
  observeEvent(input$table_group_transfo, {
    req(input$table_group_transfo)
    rv$transfoTable <- hot_to_r(input$table_group_transfo)
  })
  
  # ---------------------------------------------
  # Render design table (editable) using DT
  # ---------------------------------------------
  output$design_table <- renderDT({
    req(nrow(rv$design) > 0)
    datatable(rv$design, editable = TRUE)
  })
  
  # ---------------------------------------------
  # Export panel harmonization table to Excel
  # Aligns markers across batches and applies formatting
  # ---------------------------------------------
  output$export_panel <- downloadHandler(
    filename = function() {
      "panel_harmonization_grouped.xlsx"
    },
    content = function(file) {
      req(rv$panel_table)
      
      df <- rv$panel_table
      
      # Identify marker and fluorochrome columns
      marker_cols <- grep("_Marker$", names(df))
      fluo_cols   <- grep("_Fluo$", names(df))
      
      n_batches <- length(marker_cols)
      
      # Build list of data frames per batch
      batch_list <- lapply(seq_len(n_batches), function(i) {
        data.frame(
          Marker = df[[marker_cols[i]]],
          Fluo   = df[[fluo_cols[i]]],
          stringsAsFactors = FALSE
        )
      })
      
      # Align markers across batches
      all_markers <- unique(na.omit(unlist(lapply(batch_list, function(x) x$Marker))))
      new_df <- data.frame(matrix(nrow = length(all_markers), ncol = 0))
      
      for (i in seq_len(n_batches)) {
        tmp <- batch_list[[i]]
        aligned_marker <- rep(NA, length(all_markers))
        aligned_fluo <- rep(NA, length(all_markers))
        
        match_idx <- match(all_markers, tmp$Marker)
        valid <- !is.na(match_idx)
        
        aligned_marker[valid] <- tmp$Marker[match_idx[valid]]
        aligned_fluo[valid]   <- tmp$Fluo[match_idx[valid]]
        
        new_df <- cbind(new_df, aligned_fluo, aligned_marker)
        colnames(new_df)[(ncol(new_df)-1):ncol(new_df)] <- c(
          paste0("Batch", i, "_Fluo"),
          paste0("Batch", i, "_Marker")
        )
      }
      
      new_df$Normalize <- rep("FALSE", nrow(new_df))
      
      df <- new_df
      marker_cols <- grep("_Marker$", names(df))
      
      # Add MatchFlag to highlight markers present in multiple batches
      df$MatchFlag <- apply(df[, marker_cols], 1, function(x) sum(!is.na(x)) > 1)
      df <- df[order(-df$MatchFlag), ]
      df$MatchFlag <- NULL
      
      # Create workbook
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "Panel")
      openxlsx::writeData(wb, "Panel", df)
      
      # Apply cell locking and color styles
      lock_style <- openxlsx::createStyle(locked = TRUE)
      unlock_style <- openxlsx::createStyle(locked = FALSE)
      green_style <- openxlsx::createStyle(fgFill = "#C6EFCE")
      red_style   <- openxlsx::createStyle(fgFill = "#FFC7CE")
      
      for (col in seq_along(df)) {
        col_name <- names(df)[col]
        style <- if (grepl("_Fluo$", col_name)) lock_style else unlock_style
        openxlsx::addStyle(wb, "Panel", style, rows = 1:(nrow(df)+1), cols = col, gridExpand = TRUE)
      }
      
      # Apply conditional coloring for markers
      for (r in seq_len(nrow(df))) {
        style <- if(sum(!is.na(df[r, marker_cols])) > 1) green_style else red_style
        openxlsx::addStyle(wb, "Panel", style, rows = r+1, cols = marker_cols, gridExpand = TRUE)
      }
      
      # Data validation for Normalize column
      openxlsx::dataValidation(
        wb, "Panel", cols = ncol(df), rows = 2:(nrow(df)+1),
        type = "list", value = "TRUE,FALSE"
      )
      
      # Protect sheet
      openxlsx::protectWorksheet(wb, sheet = "Panel", protect = TRUE, password = "1234")
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  # ---------------------------------------------
  # Import panel table from Excel and align markers
  # ---------------------------------------------
  observeEvent(input$import_panel, {
    req(input$import_panel)
    
    tryCatch({
      df <- openxlsx::read.xlsx(input$import_panel$datapath)
      
      # Check minimal structure
      if (ncol(df) < 2) stop("Invalid panel file format")
      
      marker_cols <- grep("_Marker$", names(df), value = TRUE)
      fluo_cols   <- grep("_Fluo$", names(df), value = TRUE)
      
      # Map markers to rows/columns
      marker_map <- list()
      for (col in marker_cols) {
        for (r in 1:nrow(df)) {
          marker_val <- df[r, col]
          if (!is.na(marker_val) && marker_val != "") {
            if (!marker_val %in% names(marker_map)) marker_map[[marker_val]] <- list()
            marker_map[[marker_val]][[col]] <- r
          }
        }
      }
      
      # Align duplicated markers into single rows
      new_df <- df[0, ]
      processed_rows <- c()
      
      for (marker in names(marker_map)) {
        cols_with_marker <- names(marker_map[[marker]])
        if (length(cols_with_marker) > 1) {
          first_row <- df[marker_map[[marker]][[1]], ]
          new_row <- first_row
          for (col in cols_with_marker[-1]) {
            r <- marker_map[[marker]][[col]]
            new_row[[col]] <- df[r, col]
            fluo_col <- sub("_Marker$", "_Fluo", col)
            if (fluo_col %in% names(df)) new_row[[fluo_col]] <- df[r, fluo_col]
            processed_rows <- c(processed_rows, r)
          }
          processed_rows <- c(processed_rows, marker_map[[marker]][[cols_with_marker[1]]])
          new_df <- rbind(new_df, new_row)
        }
      }
      
      # Append remaining unprocessed rows
      remaining_rows <- setdiff(1:nrow(df), processed_rows)
      if (length(remaining_rows) > 0) new_df <- rbind(new_df, df[remaining_rows, ])
      
      rv$panel_table <- as.data.frame(new_df)
      showNotification("Panel imported and aligned successfully", type = "message")
      
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("Import failed:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  # =============================================
  # 1. Marker Selection & Preprocessing
  # =============================================
  observe({
    # Check that at least one batch and one FCS file exists
    if (length(rv$fs_list) >= 1 && length(rv$fs_list[[1]]) >= 1) {
      
      ff <- rv$fs_list[[1]][[1]]        # Take first FCS file
      param_data <- ff@parameters@data  # Get parameter info
      
      marker_names <- param_data$name
      marker_desc  <- param_data$desc
      
      # Create labels: "desc (name)" if description exists, else name only
      marker_labels <- ifelse(
        is.na(marker_desc) | marker_desc == "",
        marker_names,
        paste0(marker_desc, " (", marker_names, ")")
      )
      
      # Filter out FSC, SSC, and Time channels
      keep <- !grepl("^(time|fsc|ssc)", marker_names, ignore.case = TRUE)
      marker_names  <- marker_names[keep]
      marker_labels <- marker_labels[keep]
      
      # Preserve previously selected markers if defined
      selected_markers <- isolate(input$markers)
      if (is.null(selected_markers)) {
        selected_markers <- marker_names
      }
    }
  })
  
  # =============================================
  # 2. Covariate UI
  # =============================================
  output$covar_ui <- renderUI({
    req(rv$panel_table)
    
    # Get fluorochrome columns per batch
    fluo_list <- lapply(rv$fs_list, function(batch) {
      ff <- batch[[1]]  
      colnames(exprs(ff))
    })
    
    # Only keep fluorochromes common to all batches
    common_fluo <- Reduce(intersect, fluo_list)
    
    # Create selectInput for covariate
    selectInput(
      "covar",
      "There is a covariable? *In cyCombine, a covariate is a biological or experimental variable within a batch that should be preserved during batch correction.",
      choices = c("No", common_fluo)
    )
  })
  
  # =============================================
  # 3. Validate Design & Extract Data
  # =============================================
  observeEvent(input$Validate_design, {
    # Show modal during processing
    showModal(modalDialog(
      title = "Preparing data for cyCombine...",
      tags$div(style = "text-align: center;", tags$p("Please wait while your data is being processed.")),
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      # Filter panel to normalize
      panel_to_norm <- rv$panel_table[rv$panel_table$Normalize == TRUE, ]
      all_batches_df <- list()
      
      for (i in 1:input$n_batches) {
        batch_fs <- rv$fs_list[[i]]          
        batch_design <- rv$design[rv$design$Batch == i, ]
        batch_name <- get_batch_name(i)
        
        # Get fluorochromes and markers to normalize
        fluo_cols <- grep(paste0(batch_name, "_Fluo$"), names(panel_to_norm))
        fluo_names <- panel_to_norm[[paste0(batch_name, "_Fluo")]]
        markers    <- panel_to_norm[[paste0(batch_name, "_Marker")]]
        
        valid_idx <- !is.na(fluo_names)
        fluo_names <- fluo_names[valid_idx]
        markers    <- markers[valid_idx]
        
        # Extract expression matrix per FCS file
        for (j in seq_along(batch_fs)) {
          ff <- batch_fs[[j]]
          expr_mat <- exprs(ff)[, fluo_names, drop = FALSE]
          colnames(expr_mat) <- markers
          df <- as.data.frame(expr_mat)
          
          df$batch <- batch_name
          df$File  <- batch_design$File[j]
          
          # Add covariate if selected
          if (!is.null(input$covar) && input$covar != "No" && input$covar %in% colnames(expr_mat)) {
            df$covar <- expr_mat[, input$covar]
          }
          
          all_batches_df <- append(all_batches_df, list(df))
        }
      }
      
      # Combine all batches into one dataframe
      rv$unnormalized_data <- do.call(rbind, all_batches_df)
      showNotification("Data extracted and concatenated successfully.", type = "message")
      removeModal()
      
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("An error occurred:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  # =============================================
  # 4. Batch-specific Cofactor UI
  # =============================================
  output$batch_cofactor_ui <- renderUI({
    req(input$n_batches)
    
    tagList(
      lapply(1:input$n_batches, function(i) {
        batch_name <- ifelse(
          !is.null(rv$batch_names[[i]]) && rv$batch_names[[i]] != "",
          rv$batch_names[[i]],
          paste("Batch", i)
        )
        
        numericInput(
          paste0("cofactor_batch_", i),
          paste(batch_name, "cofactor"),
          value = 500,
          min = 0.1
        )
      })
    )
  })
  
  # =============================================
  # 5. Apply Transformation
  # =============================================
  observeEvent(input$apply_transform, {
    req(rv$unnormalized_data)
    
    showModal(modalDialog(
      title = "Applying Transformation",
      tags$div(style = "text-align: center;", tags$p("Markers are being transformed, please wait...")),
      footer = NULL,
      easyClose = FALSE
    ))
    
    df <- rv$unnormalized_data
    rv$isTransformed <- TRUE
    marker_cols <- setdiff(colnames(df), c("batch", "File"))
    
    # Exclude covariate column if present
    if (!is.null(df[[input$covar]])) {
      marker_cols <- setdiff(colnames(df), c("batch", "File", input$covar))
    }
    
    rv$transformed_markers <- marker_cols
    
    # Apply asinh transformation
    if (input$transform_mode == "global") {
      cofactor <- input$global_cofactor
      df[, marker_cols] <- asinh(df[, marker_cols] / cofactor)
      
      # Apply to FCS objects
      for (i in 1:input$n_batches) {
        batch_name <- get_batch_name(i)
        mapping <- get_marker_to_fluo(batch_name)
        fluo <- mapping[marker_cols]
        
        for (file in seq_along(rv$fs_list[[i]])) {
          exprs(rv$fs_list[[i]][[file]])[, fluo] <- asinh(exprs(rv$fs_list[[i]][[file]])[, fluo] / cofactor)
        }
      }
      
    } else { # Batch-specific cofactors
      for (i in 1:input$n_batches) {
        cofactor <- input[[paste0("cofactor_batch_", i)]]
        batch_name <- get_batch_name(i)
        idx <- df$batch == batch_name
        df[idx, marker_cols] <- asinh(df[idx, marker_cols] / cofactor)
        
        mapping <- get_marker_to_fluo(batch_name)
        fluo <- mapping[marker_cols]
        
        for (file in seq_along(rv$fs_list[[i]])) {
          exprs(rv$fs_list[[i]][[file]])[, fluo] <- asinh(exprs(rv$fs_list[[i]][[file]])[, fluo] / cofactor)
        }
      }
    }
    
    rv$unnormalized_data <- df
    removeModal()
    showNotification("Transformation applied to dataframe", type = "message")
  })
  
  # =============================================
  # 6. Design Validation
  # =============================================
  validate_design <- reactive({
    req(input$n_batches)
    d <- rv$design
    
    if (nrow(d) == 0) return("No files uploaded")
    
    # Each batch must have at least 1 reference
    refs_per_batch <- tapply(d$Type == "Reference", d$Batch, sum)
    if (any(is.na(refs_per_batch)) || any(refs_per_batch == 0)) return("Each batch must have at least 1 Reference")
    
    # At least 2 batches required
    if (length(unique(d$Batch)) < 2) return("At least 2 batches required")
    
    # Markers must be selected
    if (is.null(input$markers) || length(input$markers) == 0) return("Select markers")
    
    return("Design valid")
  })
  
  # =============================================
  # 7. Run cyCombine
  # =============================================
  observeEvent(input$run_cyCombine, {
    showModal(modalDialog(
      title = "Run cyCombine...",
      tags$div(style = "text-align: center;", tags$p("Please wait.")),
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      df <- rv$unnormalized_data
      marker_cols <- setdiff(colnames(df), c("batch", "File"))
      
      # Reference batch selection
      ref_batch <- input$goal_ref_batch
      if (ref_batch == "NULL") ref_batch <- NULL else ref_batch <- as.character(ref_batch)
      
      # Covariate (optional)
      covar <- if (!is.null(df[[input$covar]])) input$covar else NULL
      
      # Run batch correction
      rv$normalized_data <- batch_correct(
        df,
        covar = covar,
        ref.batch = ref_batch,
        markers = marker_cols,
        norm_method = input$norm_method, # scale or rank
        rlen = input$rlen,               # SOM iterations
        seed = 101
      )
      
      # Process normalized FCS objects
      rv$normalized_fcs <- process_normalized_fcs(rv, input)
      
      removeModal()
      showNotification("cyCombine completed successfully.", type = "message")
      
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(
        title = "Error",
        paste("An error occurred while running cyCombine:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  # =============================================
  # 1. Goal Reference Batch UI
  # =============================================
  output$goal_ref_batch_ui <- renderUI({
    req(input$n_batches)  # Ensure number of batches is defined
    
    batch_ids <- seq_len(input$n_batches)
    
    # Create choices: batch names, plus a "NULL" option
    choices <- sapply(batch_ids, function(i) get_batch_name(i))
    names(choices) <- choices
    choices <- c("NULL", choices)
    
    # Radio buttons to select goal reference batch
    radioButtons(
      inputId = "goal_ref_batch",
      label   = "Choose batch as goal reference:",
      choices = choices,
      selected = choices[2], # Default: first actual batch
      inline = TRUE
    )
  })
  
  # =============================================
  # 2. Visualization UI (Histogram)
  # =============================================
  output$visualization_ui <- renderUI({
    req(rv$normalized_data)  # Wait until data is normalized
    
    box(
      width = 12,
      tags$h3(
        "Visualization",
        style = "font-weight: bold; background-color: #fff9c4; padding: 4px; display: inline-block;"
      ),
      tags$h3("Histogram"),
      
      # Marker selection dropdown
      selectInput(
        "viz_marker_hist",
        "Select marker",
        choices = NULL
      ),
      
      # Density plot
      plotOutput("hist_plot", height = "600px", width = "600px")
    )
  })
  
  # Update the histogram marker selection dynamically
  observe({
    req(rv$normalized_data)
    
    marker_cols <- setdiff(colnames(rv$unnormalized_data), c("batch", "File"))
    
    updateSelectInput(
      session,
      "viz_marker_hist",
      choices = marker_cols
    )
  })
  
  # =============================================
  # 3. Render Histogram (Before vs After)
  # =============================================
  output$hist_plot <- renderPlot({
    req(rv$unnormalized_data, rv$normalized_data, input$viz_marker_hist)
    
    marker <- input$viz_marker_hist
    
    # Combine unnormalized and normalized data for comparison
    df <- bind_rows(
      mutate(rv$unnormalized_data, State = "Before", batch = as.factor(batch)),
      mutate(rv$normalized_data,   State = "After",  batch = as.factor(batch))
    )
    
    df$State <- factor(df$State, levels = c("Before", "After"))
    
    # Choose professional qualitative color palette
    palette <- RColorBrewer::brewer.pal(min(length(unique(df$batch)), 8), "Set2")
    
    # Density plot with faceting by "Before/After" state
    ggplot(df, aes_string(x = marker, fill = "batch", color = "batch")) +
      geom_density(alpha = 0.4, size = 1.2, adjust = 1.2) +
      facet_wrap(~State, ncol = 1, scales = "free_y") +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold", size = 13),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
        plot.margin = margin(10, 10, 10, 10)
      ) +
      labs(
        title = paste("Marker Density Before vs After cyCombine:", marker),
        x = marker,
        y = "Density",
        fill = "Batch",
        color = "Batch"
      ) +
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette)
  })
  
  # =============================================
  # 4. Process Normalized FCS Objects
  # =============================================
  process_normalized_fcs <- function(rv, input) {
    # Create new FCS objects with normalized values
    normalized_df <- rv$normalized_data
    new_fs_list <- list()
    
    marker_cols <- setdiff(colnames(normalized_df), c("batch", "File"))
    
    # Exclude covariate from markers if present
    if (!is.null(input$covar) && input$covar %in% colnames(normalized_df)) {
      marker_cols <- setdiff(colnames(normalized_df), c("batch", "File", input$covar))
    }
    
    for (i in seq_len(input$n_batches)) {
      batch_name <- get_batch_name(i)
      batch_fs <- rv$fs_list[[i]]
      batch_design <- rv$design[rv$design$Batch == i, ]
      
      # Map Marker → Fluorochrome
      fluo_names <- rv$panel_table[[paste0(batch_name, "_Fluo")]]
      markers    <- rv$panel_table[[paste0(batch_name, "_Marker")]]
      valid_idx <- !is.na(fluo_names) & !is.na(markers)
      marker_to_fluo <- setNames(fluo_names[valid_idx], markers[valid_idx])
      
      new_batch_fs <- list()
      
      for (j in seq_along(batch_fs)) {
        ff <- batch_fs[[j]]
        expr_mat <- exprs(ff)
        file_name <- batch_design$File[j]
        
        # Subset normalized data for this batch/file
        df_sub <- normalized_df[
          normalized_df$batch == batch_name &
            normalized_df$File == file_name, 
        ]
        if (nrow(df_sub) == 0) next
        
        # Replace original expression values with normalized values
        for (m in markers) {
          fluo <- marker_to_fluo[[m]]
          if (!is.null(fluo) && fluo %in% colnames(expr_mat) && m %in% colnames(df_sub)) {
            expr_mat[, fluo] <- df_sub[[m]]
          }
        }
        exprs(ff) <- expr_mat
        new_batch_fs[[j]] <- ff
      }
      
      new_fs_list[[i]] <- new_batch_fs
    }
    
    return(new_fs_list)
  }
  
  # =============================================
  # 5. Marker to Fluorochrome Mapping Helper
  # =============================================
  get_marker_to_fluo <- function(batch_name) {
    fluo <- rv$panel_table[[paste0(batch_name, "_Fluo")]]
    marker <- rv$panel_table[[paste0(batch_name, "_Marker")]]
    valid <- !is.na(fluo) & !is.na(marker)
    setNames(fluo[valid], marker[valid])
  }
  
  # =============================================
  # 6. Scatter Plot UI (Density Plot)
  # =============================================
  output$dens_plot_ui <- renderUI({
    req(rv$unnormalized_data, rv$normalized_fcs)
    
    df <- rv$unnormalized_data
    marker_cols <- setdiff(colnames(df), c("batch", "File"))
    if (!is.null(df[[input$covar]])) {
      marker_cols <- setdiff(colnames(df), c("batch", "File", input$covar))
    }
    
    fluidRow(
      box(
        width = 12,
        tags$h3("Scatter plot"),
        fluidRow(
          uiOutput("dens_file_ui"),          # File selection (assumed defined elsewhere)
          column(4, selectInput("viz_marker_x", "Marker X", choices = marker_cols)),
          column(4, selectInput("viz_marker_y", "Marker Y", choices = marker_cols))
        ),
        br(),
        plotOutput("dens_plot", height = "800px")
      )
    )
  })
  output$dens_file_ui <- renderUI({
    req(rv$fs_list, input$n_batches)
    
    n_batches <- input$n_batches
    
    fluidRow(
      lapply(1:n_batches, function(i){
        batch_name <- get_batch_name(i)
        files <- rv$design$File[rv$design$Batch == i]
        
        column(
          width = floor(12 / n_batches), 
          selectInput(
            inputId = paste0("dens_file_batch_", i),
            label = paste("File for", batch_name),
            choices = files,
            selected = files[1],
            width="50%"
          )
        )
      })
    )
  })

  # =========================
  # Render plotDens Before / After
  # =========================
  output$dens_plot <- renderPlot({
    req(input$viz_marker_x, input$viz_marker_y)
    marker_x <- input$viz_marker_x
    marker_y <- input$viz_marker_y
    n_batches <- input$n_batches
    

    all_x <- c()
    all_y <- c()
    file_indices <- vector("list", n_batches)
    
    for(i in 1:n_batches){
      batch_name <- get_batch_name(i)
      file_choice <- input[[paste0("dens_file_batch_", i)]]
      
      file_idx <- which(rv$design$File[rv$design$Batch == i] == file_choice)
      
      if(length(file_idx) == 0){
        file_indices[[i]] <- NA
        next
      } else {
        file_indices[[i]] <- file_idx[1]
      }
      
      mapping <- get_marker_to_fluo(batch_name)
      fluo_x <- mapping[[marker_x]]
      fluo_y <- mapping[[marker_y]]
      
      if(!is.null(fluo_x) && !is.null(fluo_y)){
        all_x <- c(all_x,
                   exprs(rv$fs_list[[i]][[file_idx[1]]])[, fluo_x],
                   exprs(rv$normalized_fcs[[i]][[file_idx[1]]])[, fluo_x])
        all_y <- c(all_y,
                   exprs(rv$fs_list[[i]][[file_idx[1]]])[, fluo_y],
                   exprs(rv$normalized_fcs[[i]][[file_idx[1]]])[, fluo_y])
      }
    }
    
    xlim <- range(all_x, na.rm = TRUE)
    ylim <- range(all_y, na.rm = TRUE)
    

    par(mfrow = c(2, n_batches), mar = c(3,3,2,1), pty = "s")
    

    for(i in 1:n_batches){
      file_idx <- file_indices[[i]]
      if(is.na(file_idx)) next
      
      batch_name <- get_batch_name(i)
      mapping <- get_marker_to_fluo(batch_name)
      fluo_x <- mapping[[marker_x]]
      fluo_y <- mapping[[marker_y]]
      
      if(is.null(fluo_x) || is.null(fluo_y)) next
      
      flowDensity::plotDens(
        rv$fs_list[[i]][[file_idx]],
        channels = c(fluo_x, fluo_y),
        main = paste("Before -", batch_name),
        xlab = marker_x,
        ylab = marker_y,
        xlim = xlim,
        ylim = ylim
      )
    }
    
    for(i in 1:n_batches){
      file_idx <- file_indices[[i]]
      if(is.na(file_idx)) next
      
      batch_name <- get_batch_name(i)
      mapping <- get_marker_to_fluo(batch_name)
      fluo_x <- mapping[[marker_x]]
      fluo_y <- mapping[[marker_y]]
      
      if(is.null(fluo_x) || is.null(fluo_y)) next
      
      flowDensity::plotDens(
        rv$normalized_fcs[[i]][[file_idx]],
        channels = c(fluo_x, fluo_y),
        main = paste("After -", batch_name),
        xlab = marker_x,
        ylab = marker_y,
        xlim = xlim,
        ylim = ylim
      )
    }
  })
  
  
  
  output$download_norm <- downloadHandler(
    filename = function() {
      paste0("Experiment_", input$experiment_name, "_afterCyCombine.zip")
    },
    content = function(file) {
      # Show modal to inform user
      showModal(modalDialog(
        title = "Preparing Download",
        tags$div(
          style = "text-align: center;",
          tags$p("Detransforming markers if needed and preparing ZIP file...")
        ),
        footer = NULL,
        easyClose = FALSE
      ))
      
      tmp_dir <- tempdir()
      fcs_files <- c()
      
      for (i in seq_along(rv$normalized_fcs)) {
        batch_fs <- rv$normalized_fcs[[i]]
        batch_name <- get_batch_name(i)
        
        # Récupérer mapping marker → fluo
        mapping <- get_marker_to_fluo(batch_name)
        
        for (j in seq_along(batch_fs)) {
          ff <- batch_fs[[j]]  # <-- define ff before using
          print("istransformed")
          print(rv$isTransformed)
          # Détransfo uniquement si transformation appliquée
          if (isTRUE(rv$isTransformed)) {
            marker_cols <- rv$transformed_markers
            fluo <- mapping[marker_cols]
            print(rv$transformed_markers)
            print(fluo)
            # Safety check
            fluo <- fluo[!is.na(fluo)]
            if (length(fluo) > 0) {
              if (input$transform_mode == "global") {
                cofactor <- input$global_cofactor
                exprs(ff)[, fluo] <- sinh(exprs(ff)[, fluo]) * cofactor
              } else {
                cofactor <- input[[paste0("cofactor_batch_", i)]]
                exprs(ff)[, fluo] <- sinh(exprs(ff)[, fluo]) * cofactor
              }
            }
          }
          
          # Compensation if needed
          spill <- FlowCIPHE::found.spill.CIPHE(ff)
          if (!is.null(spill) && length(spill) > 0 && !all(is.na(spill)) && spill[1] != "NULL") {
            spill_name <- spill[1]
            if (spill_name %in% names(ff@description)) {
              print(paste("Applying compensation with:", spill_name))
              ff <- flowCore::decompensate(ff, ff@description[[spill_name]])
            } else {
              print("Spill name not found in description")
            }
          } else {
            print("No valid spill matrix found")
          }
          
          # Save FCS
          name <- rv$design$File[rv$design$Batch == i][j]
          name <- gsub(".fcs", "", name)
          file_name <- paste0(batch_name, "_", name, "_norm.fcs")
          file_path <- file.path(tmp_dir, file_name)
          
          write.FCS(ff, filename = file_path)
          fcs_files <- c(fcs_files, file_path)
        }
      }
      
      # ZIP all files
      zip::zip(zipfile = file, files = fcs_files, mode = "cherry-pick")
      
      # Remove modal after done
      removeModal()
    }
  )
}


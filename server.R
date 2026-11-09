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
              truncate_max_range = FALSE,
              emptyValue = FALSE
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
    
    # ---------------------------------------------
    # Helper: normalize names for comparison
    # ---------------------------------------------
    normalize_marker <- function(x) {
      x <- as.character(x)
      x[is.na(x)] <- ""
      
      x <- trimws(tolower(x))
      x <- sub("[-_ ]+[abc]$", "", x)
      x <- gsub("[^[:alnum:]]", "", x)
      
      x
    }
    
    # ---------------------------------------------
    # Extract parameters from each batch
    # ---------------------------------------------
    for (i in seq_len(n_batches)) {
      
      req(length(rv$fs_list[[i]]) > 0)
      
      ff <- rv$fs_list[[i]][[1]]
      param_data <- ff@parameters@data
      
      fluo <- as.character(param_data$name)
      marker <- as.character(param_data$desc)
      
      # Replace empty descriptions with fluorochrome names
      missing_marker <- is.na(marker) | trimws(marker) == ""
      marker[missing_marker] <- fluo[missing_marker]
      
      # Handle mass cytometry names
      is_mass <- any(grepl(
        "^[0-9]+[A-Za-z]+_",
        marker
      ), na.rm = TRUE)
      
      if (is_mass) {
        marker <- sub(".*_", "", marker)
      }
      
      panel_list[[i]] <- data.frame(
        Fluo = fluo,
        Marker = marker,
        stringsAsFactors = FALSE
      )
    }
    
    # ---------------------------------------------
    # Batch 1 is the reference order
    # ---------------------------------------------
    reference <- panel_list[[1]]
    
    reference_marker_key <- normalize_marker(reference$Marker)
    reference_fluo_key <- normalize_marker(reference$Fluo)
    
    aligned_list <- vector("list", n_batches)
    aligned_list[[1]] <- reference
    
    # ---------------------------------------------
    # Helper: find a unique approximate match
    # ---------------------------------------------
    find_approx_match <- function(target, candidates_keys, available) {
      
      if (is.na(target) || target == "") {
        return(NA_integer_)
      }
      
      candidates <- which(
        available &
          candidates_keys != ""
      )
      
      if (length(candidates) == 0) {
        return(NA_integer_)
      }
      
      similarities <- vapply(candidates, function(j) {
        
        d <- as.numeric(utils::adist(
          target,
          candidates_keys[j]
        ))
        
        max_len <- max(
          nchar(target),
          nchar(candidates_keys[j]),
          1
        )
        
        1 - d / max_len
        
      }, numeric(1))
      
      best_score <- max(similarities)
      best <- which(similarities == best_score)
      
      threshold <- if (nchar(target) <= 4) {
        0.80
      } else {
        0.75
      }
      
      if (best_score >= threshold && length(best) == 1) {
        return(candidates[best])
      }
      
      NA_integer_
    }
    
    # ---------------------------------------------
    # Align subsequent batches to Batch 1
    # Priority:
    # 1. Exact desc
    # 2. Approximate desc
    # 3. Exact name
    # 4. Approximate name
    # ---------------------------------------------
    if (n_batches >= 2) {
      
      for (i in 2:n_batches) {
        
        current <- panel_list[[i]]
        
        current_marker_key <- normalize_marker(current$Marker)
        current_fluo_key <- normalize_marker(current$Fluo)
        
        aligned_fluo <- rep(NA_character_, nrow(reference))
        aligned_marker <- rep(NA_character_, nrow(reference))
        
        used <- rep(FALSE, nrow(current))
        
        for (r in seq_len(nrow(reference))) {
          
          matched <- FALSE
          j <- NA_integer_
          
          # ---------------------------------------
          # 1. Exact match on desc (marker)
          # ---------------------------------------
          target_marker <- reference_marker_key[r]
          
          if (!is.na(target_marker) && target_marker != "") {
            
            candidates <- which(
              !used &
                current_marker_key == target_marker &
                current_marker_key != ""
            )
            
            if (length(candidates) == 1) {
              j <- candidates[1]
            }
          }
          
          # ---------------------------------------
          # 2. Approximate match on desc
          # Only if exact desc did not match
          # ---------------------------------------
          if (is.na(j)) {
            
            j <- find_approx_match(
              target_marker,
              current_marker_key,
              !used
            )
          }
          
          # If desc matched, store and continue
          if (!is.na(j)) {
            
            aligned_fluo[r] <- current$Fluo[j]
            aligned_marker[r] <- current$Marker[j]
            used[j] <- TRUE
            
            next
          }
          
          # ---------------------------------------
          # 3. Exact match on name (fluorochrome)
          # Fallback if desc did not match
          # ---------------------------------------
          target_fluo <- reference_fluo_key[r]
          
          if (!is.na(target_fluo) && target_fluo != "") {
            
            candidates <- which(
              !used &
                current_fluo_key == target_fluo &
                current_fluo_key != ""
            )
            
            if (length(candidates) == 1) {
              j <- candidates[1]
            }
          }
          
          # ---------------------------------------
          # 4. Approximate match on name
          # ---------------------------------------
          if (is.na(j)) {
            
            j <- find_approx_match(
              target_fluo,
              current_fluo_key,
              !used
            )
          }
          
          if (!is.na(j)) {
            
            aligned_fluo[r] <- current$Fluo[j]
            aligned_marker[r] <- current$Marker[j]
            used[j] <- TRUE
          }
        }
        
        # -----------------------------------------
        # Append unmatched parameters at the end
        # -----------------------------------------
        unmatched <- which(!used)
        
        aligned <- data.frame(
          Fluo = c(
            aligned_fluo,
            current$Fluo[unmatched]
          ),
          Marker = c(
            aligned_marker,
            current$Marker[unmatched]
          ),
          stringsAsFactors = FALSE
        )
        
        aligned_list[[i]] <- aligned
      }
    }
    
    # ---------------------------------------------
    # Build final table
    # Batch 1 order is preserved
    # ---------------------------------------------
    max_len <- max(vapply(
      aligned_list,
      nrow,
      integer(1)
    ))
    
    panel_table <- data.frame(
      matrix(nrow = max_len, ncol = 0)
    )
    
    for (i in seq_len(n_batches)) {
      
      df <- aligned_list[[i]]
      
      if (nrow(df) < max_len) {
        
        extra <- max_len - nrow(df)
        
        df <- rbind(
          df,
          data.frame(
            Fluo = rep(NA_character_, extra),
            Marker = rep(NA_character_, extra),
            stringsAsFactors = FALSE
          )
        )
      }
      
      batch_name <- get_batch_name(i)
      
      colnames(df) <- c(
        paste0(batch_name, "_Fluo"),
        paste0(batch_name, "_Marker")
      )
      
      panel_table <- cbind(panel_table, df)
    }
    
    
    # Conserver les sélections existantes lors de la reconstruction
    old_panel <- isolate(rv$panel_table)
    
    # Par défaut, aucune sélection si c'est la première création
    panel_table$Normalize <- FALSE
    
    if (!is.null(old_panel) &&
        "Normalize" %in% names(old_panel)) {
      
      ref_fluo_col <- paste0(get_batch_name(1), "_Fluo")
      ref_marker_col <- paste0(get_batch_name(1), "_Marker")
      
      old_cols_ok <- all(
        c(ref_fluo_col, ref_marker_col) %in% names(old_panel)
      )
      
      new_cols_ok <- all(
        c(ref_fluo_col, ref_marker_col) %in% names(panel_table)
      )
      
      if (old_cols_ok && new_cols_ok) {
        
        # Construire une clé uniquement pour les lignes valides
        old_valid <- (
          !is.na(old_panel[[ref_fluo_col]]) &
            !is.na(old_panel[[ref_marker_col]]) &
            trimws(as.character(old_panel[[ref_fluo_col]])) != "" &
            trimws(as.character(old_panel[[ref_marker_col]])) != ""
        )
        
        new_valid <- (
          !is.na(panel_table[[ref_fluo_col]]) &
            !is.na(panel_table[[ref_marker_col]]) &
            trimws(as.character(panel_table[[ref_fluo_col]])) != "" &
            trimws(as.character(panel_table[[ref_marker_col]])) != ""
        )
        
        old_key <- paste(
          as.character(old_panel[[ref_fluo_col]]),
          as.character(old_panel[[ref_marker_col]]),
          sep = "|||"
        )
        
        new_key <- paste(
          as.character(panel_table[[ref_fluo_col]]),
          as.character(panel_table[[ref_marker_col]]),
          sep = "|||"
        )
        
        old_normalize <- tolower(
          trimws(as.character(old_panel$Normalize))
        )
        
        selected_keys <- old_key[
          old_valid &
            !is.na(old_normalize) &
            old_normalize %in% c("true", "t", "1", "yes")
        ]
        
        panel_table$Normalize <- (
          new_valid & new_key %in% selected_keys
        )
      }
    }
    
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
      
      # 
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
  
  # =============================================
  # Validate Design & Extract Data
  # Noms de référence = marqueurs du Batch 1
  # Les fichiers FCS originaux restent inchangés
  # =============================================
  observeEvent(input$Validate_design, {
    
    showModal(modalDialog(
      title = "Preparing data for cyCombine...",
      tags$p("Please wait while your data is being processed."),
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      
      req(rv$panel_table)
      req(input$n_batches)
      
      panel <- rv$panel_table
      all_batches_df <- list()
      
      # Les lignes cochées pour normalisation
      normalize_rows <- which(
        !is.na(panel$Normalize) & panel$Normalize == TRUE
      )
      
      if (length(normalize_rows) == 0) {
        stop("Please select at least one marker to normalize.")
      }
      
      # Noms de référence du Batch 1
      ref_col <- paste0(get_batch_name(1), "_Marker")
      if (!ref_col %in% names(panel)) {
        stop("Reference marker column for Batch 1 not found.")
      }
      
      reference_markers <- as.character(panel[[ref_col]][normalize_rows])
      
      if (anyNA(reference_markers) ||
          any(trimws(reference_markers) == "")) {
        stop("Some selected markers have no reference name in Batch 1.")
      }
      
      if (anyDuplicated(reference_markers)) {
        stop("Selected reference marker names in Batch 1 must be unique.")
      }
      
      for (i in seq_len(input$n_batches)) {
        
        batch_name <- get_batch_name(i)
        batch_fs <- rv$fs_list[[i]]
        batch_design <- rv$design[rv$design$Batch == i, ]
        
        fluo_col <- paste0(batch_name, "_Fluo")
        marker_col <- paste0(batch_name, "_Marker")
        
        if (!fluo_col %in% names(panel) ||
            !marker_col %in% names(panel)) {
          stop(paste("Panel columns not found for", batch_name))
        }
        
        # Correspondance par ligne de la table d'harmonisation.
        # On ne fait PAS de jointure sur les noms de marqueurs.
        fluo_names <- as.character(panel[[fluo_col]][normalize_rows])
        
        if (anyNA(fluo_names) ||
            any(trimws(fluo_names) == "")) {
          stop(paste(
            "A selected marker has no corresponding channel in",
            batch_name
          ))
        }
        
        for (j in seq_along(batch_fs)) {
          
          ff <- batch_fs[[j]]
          mat <- flowCore::exprs(ff)
          
          if (anyDuplicated(fluo_names)) {
            stop(paste("Duplicate channels selected in", batch_name))
          }
          
          missing_channels <- setdiff(fluo_names, colnames(mat))
          
          if (length(missing_channels) > 0) {
            stop(paste(
              "Channels missing in", batch_name, ":",
              paste(missing_channels, collapse = ", ")
            ))
          }
          
          expr_mat <- mat[, fluo_names, drop = FALSE]
          
          # Même nom backend dans tous les batches :
          # le nom de référence du Batch 1.
          colnames(expr_mat) <- reference_markers
          
          df <- as.data.frame(expr_mat, check.names = FALSE)
          df$batch <- batch_name
          df$File <- batch_design$File[j]
          
          all_batches_df[[length(all_batches_df) + 1L]] <- df
        }
      }
      
      rv$unnormalized_data <- do.call(rbind, all_batches_df)
      
      showNotification(
        "Data extracted with Batch 1 marker names as backend reference.",
        type = "message"
      )
      
      removeModal()
      
    }, error = function(e) {
      
      removeModal()
      
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
    req(rv$panel_table)
    req(rv$fs_list)
    
    showModal(modalDialog(
      title = "Applying Transformation",
      tags$div(
        style = "text-align: center;",
        tags$p("Markers are being transformed, please wait...")
      ),
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      
      df <- rv$unnormalized_data
      
      # Marqueurs communs utilisés en backend
      marker_cols <- setdiff(
        colnames(df),
        c("batch", "File")
      )
      
      # Exclure la covariable si elle existe réellement
      covar <- input$covar
      
      if (!is.null(covar) &&
          length(covar) == 1 &&
          !is.na(covar) &&
          nzchar(covar) &&
          covar %in% colnames(df)) {
        
        marker_cols <- setdiff(marker_cols, covar)
      }
      
      validate(
        need(length(marker_cols) > 0,
             "Aucun marqueur à transformer.")
      )
      
      rv$transformed_markers <- marker_cols
      
      # Noms de référence du Batch 1
      reference_batch <- get_batch_name(1)
      reference_col <- paste0(reference_batch, "_Marker")
      
      validate(
        need(reference_col %in% colnames(rv$panel_table),
             "Colonne des marqueurs du Batch 1 introuvable.")
      )
      
      reference_markers <- trimws(
        as.character(rv$panel_table[[reference_col]])
      )
      
      panel <- rv$panel_table
      
      # Appliquer la transformation au dataframe et aux objets FCS en mémoire
      for (i in seq_len(input$n_batches)) {
        
        batch_name <- get_batch_name(i)
        fluo_col <- paste0(batch_name, "_Fluo")
        
        validate(
          need(fluo_col %in% colnames(panel),
               paste("Colonne des canaux introuvable :", fluo_col))
        )
        
        fluo_names <- trimws(as.character(panel[[fluo_col]]))
        
        valid_idx <- (
          !is.na(reference_markers) &
            reference_markers != "" &
            !is.na(fluo_names) &
            fluo_names != ""
        )
        
        marker_to_fluo <- stats::setNames(
          fluo_names[valid_idx],
          reference_markers[valid_idx]
        )
        
        # Ne garder que les marqueurs demandés présents dans ce batch
        current_markers <- intersect(
          marker_cols,
          names(marker_to_fluo)
        )
        
        if (length(current_markers) == 0) {
          next
        }
        
        current_fluo <- unname(
          marker_to_fluo[current_markers]
        )
        
        # Vérifier les canaux dans tous les fichiers du batch
        for (file in seq_along(rv$fs_list[[i]])) {
          
          fcs_obj <- rv$fs_list[[i]][[file]]
          mat <- flowCore::exprs(fcs_obj)
          
          missing_fluo <- setdiff(current_fluo, colnames(mat))
          
          if (length(missing_fluo) > 0) {
            stop(
              "Canal(s) introuvable(s) dans le batch ",
              batch_name, ", fichier ", file, " : ",
              paste(missing_fluo, collapse = ", ")
            )
          }
          
          # Déterminer les lignes du dataframe pour ce fichier
          design_rows <- which(rv$design$Batch == i)
          
          if (length(design_rows) < file) {
            stop(
              "Impossible d'associer le fichier ", file,
              " au batch ", batch_name, " dans rv$design."
            )
          }
          
          design_row <- design_rows[file]
          file_name <- as.character(rv$design$File[design_row])
          batch_label <- as.character(rv$design$BatchName[design_row])
          
          df_rows <- (
            as.character(df$batch) == batch_label &
              as.character(df$File) == file_name
          )
          
          df_rows[is.na(df_rows)] <- FALSE
          
          # Extraire le cofacteur
          if (input$transform_mode == "global") {
            
            cofactor <- input$global_cofactor
            
          } else {
            
            cofactor <- input[[paste0("cofactor_batch_", i)]]
          }
          
          validate(
            need(
              length(cofactor) == 1 &&
                is.finite(cofactor) &&
                cofactor > 0,
              paste("Cofacteur invalide pour le batch", batch_name)
            )
          )
          
          # Transformer les valeurs du dataframe
          if (any(df_rows)) {
            
            if (sum(df_rows) != nrow(mat)) {
              stop(
                "Le nombre d'événements ne correspond pas pour ",
                file_name, " : dataframe = ", sum(df_rows),
                ", FCS = ", nrow(mat), "."
              )
            }
            
            df[df_rows, current_markers] <- lapply(
              df[df_rows, current_markers, drop = FALSE],
              function(x) asinh(as.matrix(x) / cofactor)
            )
          }
          
          # Transformer les canaux correspondants dans l'objet FCS
          mat[, current_fluo] <- asinh(
            mat[, current_fluo, drop = FALSE] / cofactor
          )
          
          flowCore::exprs(fcs_obj) <- mat
          rv$fs_list[[i]][[file]] <- fcs_obj
        }
      }
      
      rv$unnormalized_data <- df
      rv$isTransformed <- TRUE
      
      removeModal()
      
      showNotification(
        "Transformation appliquée au dataframe et aux objets FCS en mémoire.",
        type = "message"
      )
      
    }, error = function(e) {
      
      removeModal()
      
      showNotification(
        paste("Erreur pendant la transformation :", conditionMessage(e)),
        type = "error",
        duration = NULL
      )
      
    })
    
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
      
      print(df)
      print(input$goal_ref_batch)
      
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
      print("ok")
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
      
      df <- bind_rows(
        mutate(rv$unnormalized_data, State = "Before", batch = as.factor(batch)),
        mutate(rv$normalized_data,   State = "After",  batch = as.factor(batch))
      )
      
      df$State <- factor(df$State, levels = c("Before", "After"))
      
      # palette safe
      nbatch <- length(unique(df$batch))
      palette <- RColorBrewer::brewer.pal(min(max(nbatch, 3), 8), "Set2")
      
      ggplot(df, aes(x = .data[[marker]], fill = batch, color = batch)) +
        geom_density(alpha = 0.4, linewidth = 1.1, adjust = 1.2) +
        facet_wrap(~State, ncol = 1, scales = "free_y") +
        theme_minimal(base_size = 14) +
        theme(
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank(),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          legend.title = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey95", color = NA),
          strip.text = element_text(face = "bold")
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
    
    # Vérifications de base
    if (is.null(rv$panel_table) ||
        is.null(rv$normalized_data) ||
        is.null(rv$fs_list)) {
      stop("Les données normalisées, le panel ou les fichiers FCS sont manquants.")
    }
    
    normalized_df <- rv$normalized_data
    panel <- rv$panel_table
    
    if (!all(c("batch", "File") %in% colnames(normalized_df))) {
      stop("Les colonnes 'batch' et 'File' sont absentes des données normalisées.")
    }
    
    # Les noms du Batch 1 servent uniquement de noms internes communs
    reference_batch_name <- get_batch_name(1)
    reference_col <- paste0(reference_batch_name, "_Marker")
    
    if (!reference_col %in% colnames(panel)) {
      stop(
        "Colonne de marqueurs de référence introuvable : ",
        reference_col
      )
    }
    
    reference_markers <- trimws(as.character(panel[[reference_col]]))
    
    # Préparer une copie des objets FCS en mémoire.
    # Les fichiers source sur disque ne sont jamais modifiés.
    normalized_fcs <- rv$fs_list
    
    for (i in seq_len(input$n_batches)) {
      
      batch_name <- get_batch_name(i)
      fluo_col <- paste0(batch_name, "_Fluo")
      
      if (!fluo_col %in% colnames(panel)) {
        stop("Colonne de canaux introuvable : ", fluo_col)
      }
      
      batch_fs <- rv$fs_list[[i]]
      fluo_names <- trimws(as.character(panel[[fluo_col]]))
      
      # Correspondance par ligne du panel :
      # nom de référence du Batch 1 -> canal réel du batch courant.
      valid_idx <- (
        !is.na(reference_markers) &
          reference_markers != "" &
          !is.na(fluo_names) &
          fluo_names != ""
      )
      
      marker_to_fluo <- stats::setNames(
        fluo_names[valid_idx],
        reference_markers[valid_idx]
      )
      
      if (anyDuplicated(names(marker_to_fluo))) {
        stop(
          "Les noms de marqueurs de référence du Batch 1 ne sont pas uniques."
        )
      }
      
      if (anyDuplicated(unname(marker_to_fluo))) {
        stop(
          "Plusieurs marqueurs sélectionnés correspondent au même canal dans ",
          batch_name, "."
        )
      }
      
      # Vérifier que les canaux existent réellement dans les FCS du batch.
      available_channels <- as.character(
        batch_fs[[1]]@parameters@data$name
      )
      
      missing_channels <- setdiff(
        unname(marker_to_fluo),
        available_channels
      )
      
      if (length(missing_channels) > 0) {
        stop(
          "Canaux absents du batch ", batch_name, " : ",
          paste(missing_channels, collapse = ", ")
        )
      }
      
      # Colonnes normalisées communes, nommées selon le Batch 1.
      markers_to_update <- intersect(
        names(marker_to_fluo),
        colnames(normalized_df)
      )
      
      if (length(markers_to_update) == 0) {
        next
      }
      
      for (j in seq_along(batch_fs)) {
        
        # Retrouver les événements correspondant au fichier.
        file_name <- as.character(
          rv$design$File[
            rv$design$Batch == i
          ][j]
        )
        
        batch_label <- as.character(
          rv$design$BatchName[
            rv$design$Batch == i
          ][j]
        )
        
        rows <- (
          as.character(normalized_df$batch) == batch_label &
            as.character(normalized_df$File) == file_name
        )
        
        rows[is.na(rows)] <- FALSE
        df_sub <- normalized_df[rows, , drop = FALSE]
        
        if (nrow(df_sub) == 0) {
          warning(
            "Aucune donnée normalisée trouvée pour ",
            batch_name, " / ", file_name,
            ". Le fichier est conservé sans modification."
          )
          next
        }
        
        fcs_obj <- batch_fs[[j]]
        expr_mat <- flowCore::exprs(fcs_obj)
        
        if (nrow(df_sub) != nrow(expr_mat)) {
          stop(
            "Le nombre d'événements ne correspond pas pour ",
            file_name, " (batch ", batch_name, "). ",
            "FCS : ", nrow(expr_mat),
            " ; données normalisées : ", nrow(df_sub), "."
          )
        }
        
        # Réinjecter chaque marqueur normalisé dans son canal d'origine.
        for (marker in markers_to_update) {
          
          fluo <- unname(marker_to_fluo[marker])
          
          if (!fluo %in% colnames(expr_mat)) {
            stop(
              "Le canal '", fluo, "' est absent de ",
              file_name, "."
            )
          }
          
          if (!marker %in% colnames(df_sub)) {
            next
          }
          
          values <- df_sub[[marker]]
          
          if (length(values) != nrow(expr_mat)) {
            stop(
              "Le nombre de valeurs ne correspond pas pour le marqueur '",
              marker, "' dans ", file_name, "."
            )
          }
          
          expr_mat[, fluo] <- as.numeric(values)
        }
        
        # Modifier uniquement l'objet FCS en mémoire.
        flowCore::exprs(fcs_obj) <- expr_mat
        batch_fs[[j]] <- fcs_obj
      }
      
      normalized_fcs[[i]] <- batch_fs
    }
    
    return(normalized_fcs)
  }
  
  # =============================================
  # 5. Marker to Fluorochrome Mapping Helper
  # =============================================
  get_marker_to_fluo <- function(batch_name) {
    
    reference_batch <- get_batch_name(1)
    
    reference_col <- paste0(reference_batch, "_Marker")
    fluo_col <- paste0(batch_name, "_Fluo")
    
    if (!all(c(reference_col, fluo_col) %in% names(rv$panel_table))) {
      return(stats::setNames(character(0), character(0)))
    }
    
    reference_markers <- trimws(
      as.character(rv$panel_table[[reference_col]])
    )
    
    fluo_names <- trimws(
      as.character(rv$panel_table[[fluo_col]])
    )
    
    valid <- (
      !is.na(reference_markers) &
        reference_markers != "" &
        !is.na(fluo_names) &
        fluo_names != ""
    )
    
    mapping <- stats::setNames(
      fluo_names[valid],
      reference_markers[valid]
    )
    
    mapping
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
    
    req(
      input$viz_marker_x,
      input$viz_marker_y,
      rv$fs_list,
      rv$normalized_fcs,
      rv$panel_table
    )
    
    marker_x <- input$viz_marker_x
    marker_y <- input$viz_marker_y
    n_batches <- min(
      input$n_batches,
      length(rv$fs_list),
      length(rv$normalized_fcs)
    )
    
    validate(
      need(n_batches > 0, "Aucun batch disponible."),
      need(marker_x != marker_y, "Choisis deux marqueurs différents.")
    )
    
    # Vérifier les correspondances et déterminer les limites communes
    file_indices <- vector("list", n_batches)
    mappings <- vector("list", n_batches)
    
    all_x <- numeric(0)
    all_y <- numeric(0)
    
    for (i in seq_len(n_batches)) {
      
      batch_name <- get_batch_name(i)
      
      file_choice <- input[[paste0("dens_file_batch_", i)]]
      batch_files <- rv$design$File[rv$design$Batch == i]
      
      file_idx <- match(file_choice, batch_files)
      
      if (is.na(file_idx)) {
        file_indices[[i]] <- NA_integer_
        next
      }
      
      file_indices[[i]] <- file_idx
      
      mapping <- get_marker_to_fluo(batch_name)
      mappings[[i]] <- mapping
      
      fluo_x <- unname(mapping[marker_x])
      fluo_y <- unname(mapping[marker_y])
      
      if (
        length(fluo_x) != 1 || is.na(fluo_x) ||
        length(fluo_y) != 1 || is.na(fluo_y) ||
        !nzchar(fluo_x) || !nzchar(fluo_y)
      ) {
        next
      }
      
      # Avant
      mat_before <- flowCore::exprs(
        rv$fs_list[[i]][[file_idx]]
      )
      
      if (!all(c(fluo_x, fluo_y) %in% colnames(mat_before))) {
        next
      }
      
      # Après
      mat_after <- flowCore::exprs(
        rv$normalized_fcs[[i]][[file_idx]]
      )
      
      if (!all(c(fluo_x, fluo_y) %in% colnames(mat_after))) {
        next
      }
      
      all_x <- c(
        all_x,
        mat_before[, fluo_x],
        mat_after[, fluo_x]
      )
      
      all_y <- c(
        all_y,
        mat_before[, fluo_y],
        mat_after[, fluo_y]
      )
    }
    
    validate(
      need(length(all_x) > 0 && length(all_y) > 0,
           "Aucune donnée disponible pour ces marqueurs. Vérifie leur correspondance dans le panel.")
    )
    
    xlim <- range(all_x[is.finite(all_x)], na.rm = TRUE)
    ylim <- range(all_y[is.finite(all_y)], na.rm = TRUE)
    
    validate(
      need(all(is.finite(xlim)) && diff(xlim) > 0,
           "Les valeurs du marqueur X ne permettent pas de tracer le graphique."),
      need(all(is.finite(ylim)) && diff(ylim) > 0,
           "Les valeurs du marqueur Y ne permettent pas de tracer le graphique.")
    )
    
    # Disposition : première ligne = avant, seconde ligne = après
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    
    par(
      mfrow = c(2, n_batches),
      mar = c(4, 4, 3, 1),
      pty = "s"
    )
    
    # Graphiques avant normalisation
    for (i in seq_len(n_batches)) {
      
      file_idx <- file_indices[[i]]
      batch_name <- get_batch_name(i)
      
      if (is.na(file_idx)) {
        plot.new()
        title(main = paste("Before -", batch_name))
        text(0.5, 0.5, "Aucun fichier sélectionné")
        next
      }
      
      mapping <- mappings[[i]]
      fluo_x <- unname(mapping[marker_x])
      fluo_y <- unname(mapping[marker_y])
      
      if (
        length(fluo_x) != 1 || is.na(fluo_x) ||
        length(fluo_y) != 1 || is.na(fluo_y) ||
        !nzchar(fluo_x) || !nzchar(fluo_y)
      ) {
        plot.new()
        title(main = paste("Before -", batch_name))
        text(0.5, 0.5, "Marqueur non disponible")
        next
      }
      
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
    
    # Graphiques après normalisation
    for (i in seq_len(n_batches)) {
      
      file_idx <- file_indices[[i]]
      batch_name <- get_batch_name(i)
      
      if (is.na(file_idx)) {
        plot.new()
        title(main = paste("After -", batch_name))
        text(0.5, 0.5, "Aucun fichier sélectionné")
        next
      }
      
      mapping <- mappings[[i]]
      fluo_x <- unname(mapping[marker_x])
      fluo_y <- unname(mapping[marker_y])
      
      if (
        length(fluo_x) != 1 || is.na(fluo_x) ||
        length(fluo_y) != 1 || is.na(fluo_y) ||
        !nzchar(fluo_x) || !nzchar(fluo_y)
      ) {
        plot.new()
        title(main = paste("After -", batch_name))
        text(0.5, 0.5, "Marqueur non disponible")
        next
      }
      
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
        
        # R??cup??rer mapping marker ??? fluo
        mapping <- get_marker_to_fluo(batch_name)
        
        for (j in seq_along(batch_fs)) {
          ff <- batch_fs[[j]] 
          print("istransformed")
          print(rv$isTransformed)
          #
          if (isTRUE(rv$isTransformed)) {
            marker_cols <- rv$transformed_markers
            fluo <- mapping[marker_cols]
            print(rv$transformed_markers)
            print(fluo)
            
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


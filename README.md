# cyCombine Shiny Implementation

## Overview

This repository contains a Shiny implementation of the existing **cyCombine** tool for multi-batch cytometry data integration and normalization (PMID: 35361793).

Original cyCombine repository: https://github.com/biosurf/cyCombine

This Shiny application provides a user-friendly interface to:

- Upload multiple batches of FCS files.
- Automatically detect and apply compensation matrices when available.
- Harmonize panels across batches.
- Select markers for normalization.
- Apply an arcsinh transformation with a global or batch-specific cofactor.
- Perform batch correction using cyCombine.
- Visualize marker distributions before and after normalization.
- Compare cell distributions using density plots.
- Export normalized FCS files.

For detailed instructions, please refer to the accompanying PDF tutorial.

---

## 1. Requirements

Before running the application, make sure you have installed:

- **R** (https://cran.r-project.org/)
- **RStudio** (recommended): https://posit.co/download/rstudio-desktop/

An internet connection is required to install the R packages.

---

## 2. Installation of required R packages

The application requires packages from CRAN, Bioconductor and GitHub.

### 2.1. Install CRAN packages

Run the following commands in the R console:

```r
install.packages(c(
  "shiny",
  "shinydashboard",
  "shinyjs",
  "shinycssloaders",
  "rhandsontable",
  "DT",
  "dplyr",
  "ggplot2",
  "openxlsx",
  "flowDensity",
  "RColorBrewer",
  "zip",
  "remotes"
))
```

### 2.2. Install Bioconductor packages

The application uses `flowCore` to read, manipulate and export FCS files.

The original cyCombine package also requires Bioconductor dependencies for cytometry data correction.

Install them using:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "flowCore",
  "Biobase",
  "sva"
))
```

### 2.3. Install cyCombine from GitHub

Install the original cyCombine package:

```r
remotes::install_github("biosurf/cyCombine")
```

### 2.4. Install FlowCIPHE from GitHub

The application uses FlowCIPHE to automatically detect compensation matrices stored in FCS files.

Install it using:

```r
remotes::install_github("cipheLab/FlowCIPHE")
```

---

## 3. Check package installation

To check that the main packages are correctly installed, run:

```r
required_packages <- c(
  "shiny",
  "shinydashboard",
  "shinyjs",
  "shinycssloaders",
  "rhandsontable",
  "DT",
  "flowCore",
  "flowDensity",
  "cyCombine",
  "FlowCIPHE",
  "dplyr",
  "ggplot2",
  "openxlsx",
  "RColorBrewer",
  "zip"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) == 0) {
  message("All required packages are installed.")
} else {
  message(
    "Missing packages: ",
    paste(missing_packages, collapse = ", ")
  )
}
```

---

## 4. Run the Shiny application

1. Download or clone this GitHub repository.
2. Open the project folder in RStudio.
3. Make sure all required packages are installed.
4. Open `ui.R` and `server.R`.
5. Launch the application by running:

```r
shiny::runApp()
```

Alternatively, if the files are in a specific directory:

```r
shiny::runApp("path/to/cyCombine")
```

---

## 5. Application workflow

The application is organized into two main tabs.

### Data

- Upload FCS files for each batch.
- Define batch names.
- Harmonize marker panels.
- Export and import the panel harmonization table.
- Select markers to normalize.
- Choose a biological or experimental covariate.
- Apply an arcsinh transformation.

### Run cyCombine

- Select the reference batch.
- Choose the normalization method (`scale` or `rank`).
- Set the number of SOM iterations.
- Run cyCombine for batch correction.
- Visualize marker density distributions before and after normalization.
- Compare two markers using density plots.
- Download normalized FCS files in a ZIP archive.

---

## 6. References

**cyCombine**

Dam, S. H. et al. (2022). Robust integration of single-cell cytometry datasets within and across technologies. *Nature Communications*.

GitHub: https://github.com/biosurf/cyCombine

**FlowCIPHE**

GitHub: https://github.com/cipheLab/FlowCIPHE

---

## 7. Tutorial

For a detailed description of the application and its different functionalities, please refer to the accompanying PDF tutorial:

`cyCombine app tutorial.pdf`
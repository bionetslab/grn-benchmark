# install_dependencies.R

# Helper function to install CRAN and Bioconductor packages
install_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      message(paste("Installing package:", pkg))
      install.packages(pkg, dependencies = TRUE)
    } else {
      message(paste("Package already installed:", pkg))
    }
  }
}

# Load or install BiocManager for installing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  library(BiocManager)
}

# Bioconductor packages to install
bioconductor_packages <- c(
  "pcaMethods", "multtest", "GEOquery", "affy", 
  "genefilter", "GOstats", "ath1121501.db"
)

# CRAN packages to install
cran_packages <- c("spatstat", "igraph", "cluster", "DiffCorr")

# Install Bioconductor packages
BiocManager::install(bioconductor_packages, force = TRUE)

# Install CRAN packages
install_packages(cran_packages)

# Verify installations
message("Verifying package installations...")
required_packages <- c(
  "pcaMethods", "multtest", "GEOquery", "affy", 
  "genefilter", "GOstats", "ath1121501.db", 
  "spatstat", "igraph", "cluster", "DiffCorr"
)

# Loop through and verify
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    stop(paste("Package not installed:", pkg))
  }
}

message("Installation done")

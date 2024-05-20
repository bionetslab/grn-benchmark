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
bioconductor_packages <- c("WGCNA", "impute", "preprocessCore", "minet", "org.Mm.eg.db", "clusterProfiler", "org.Dm.eg.db", "AnnotationDbi", "enrichplot", "pbapply")

# CRAN packages to install
cran_packages <- c("dplyr", "readr", "optparse", "tidyr", "ggplot2")

# Install Bioconductor packages
BiocManager::install(bioconductor_packages)

# Install CRAN packages
install_packages(cran_packages)

# Install DGCA from archive
install.packages("https://cran.r-project.org/src/contrib/Archive/DGCA/DGCA_1.0.2.tar.gz", repos = NULL, type = "source")

# Verify installations
message("Verifying package installations...")
required_packages <- c("dplyr", "readr", "optparse", "DGCA",
                       "WGCNA", "impute", "preprocessCore", "minet", "org.Mm.eg.db",
                       "clusterProfiler", "org.Dm.eg.db", "enrichplot", "AnnotationDbi",
                       "ggplot2", "tidyr")

# Loop through and verify
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    stop(paste("Package not installed:", pkg))
  }
}

message("All packages installed successfully.")

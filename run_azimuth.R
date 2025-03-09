
library(Seurat)
library(Azimuth)
library(SeuratData)
library(reticulate)
reticulate::py_install("anndata")
reticulate::py_install("numpy")
numpy <- reticulate::import("numpy")
anndata <- import("anndata")

available_data<-AvailableData()
available_data

InstallData('kidneyref')

# Set paths
input_file <- "data/expr.h5ad"
output_file <- "data/annotated_expr.h5ad"

# Run Azimuth for cell type annotation
seurat_obj <- Azimuth::RunAzimuth(query = input_file, reference = "kidneyref")

# Inspect the first few rows of the metadata to check for 'predicted.id'
head(seurat_obj@meta.data)

# Check the available columns in `meta.data` to confirm the predicted cell type columns
colnames(seurat_obj@meta.data)

# Extract predictions for `l1`, `l2`, and `l3` from `meta.data`
predicted_l1 <- seurat_obj@meta.data$predicted.annotation.l1
predicted_l2 <- seurat_obj@meta.data$predicted.annotation.l2
predicted_l3 <- seurat_obj@meta.data$predicted.annotation.l3

# Check the first few predictions for each level
head(predicted_l1)
head(predicted_l2)
head(predicted_l3)

# Extract the counts matrix from the 'RNA' assay (stored in layers)
counts_matrix <- seurat_obj@assays$RNA@layers$counts

# Convert the counts matrix to a numpy array (Python format)
counts_matrix_python <- numpy$array(as.matrix(counts_matrix))

# Create an AnnData object using the counts matrix
adata <- anndata$AnnData(X = counts_matrix_python)

# Trim predictions to match the number of cells
predicted_l1 <- predicted_l1[1:nrow(adata$obs)]
predicted_l2 <- predicted_l2[1:nrow(adata$obs)]
predicted_l3 <- predicted_l3[1:nrow(adata$obs)]

# Add predictions to AnnData object
adata$obs$predicted_l1 <- predicted_l1
adata$obs$predicted_l2 <- predicted_l2
adata$obs$predicted_l3 <- predicted_l3

# Save the AnnData object to a .h5ad file
adata$write("data/annotated_expr_l2.h5ad")

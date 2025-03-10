
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

# Extract predictions and other metadata
predicted_l1 <- seurat_obj@meta.data$predicted.annotation.l1
predicted_l2 <- seurat_obj@meta.data$predicted.annotation.l2
predicted_l3 <- seurat_obj@meta.data$predicted.annotation.l3
predicted_l1_score <- seurat_obj@meta.data$predicted.annotation.l1.score
predicted_l2_score <- seurat_obj@meta.data$predicted.annotation.l2.score
predicted_l3_score <- seurat_obj@meta.data$predicted.annotation.l3.score
mapping_score <- seurat_obj@meta.data$mapping.score

# Extract the counts matrix from the 'RNA' assay (stored in layers)
counts_matrix <- seurat_obj@assays$RNA@layers$counts

# Convert the counts matrix to a numpy array (Python format)
counts_matrix_python <- numpy$array(as.matrix(counts_matrix))

# Create an AnnData object using the counts matrix
adata <- anndata$AnnData(X = counts_matrix_python)

# Trim predictions to match the number of cells in AnnData object
num_cells_adata <- nrow(adata$obs)  # Get the number of cells in the AnnData object
predicted_l1 <- predicted_l1[1:num_cells_adata]  # Trim predictions to match the number of cells
predicted_l2 <- predicted_l2[1:num_cells_adata]
predicted_l3 <- predicted_l3[1:num_cells_adata]
predicted_l1_score <- predicted_l1_score[1:num_cells_adata]
predicted_l2_score <- predicted_l2_score[1:num_cells_adata]
predicted_l3_score <- predicted_l3_score[1:num_cells_adata]
mapping_score <- mapping_score[1:num_cells_adata]

# Add trimmed predictions to the AnnData object
adata$obs$predicted_l1 <- predicted_l1
adata$obs$predicted_l2 <- predicted_l2
adata$obs$predicted_l3 <- predicted_l3
adata$obs$predicted_l1_score <- predicted_l1_score
adata$obs$predicted_l2_score <- predicted_l2_score
adata$obs$predicted_l3_score <- predicted_l3_score
adata$obs$mapping_score <- mapping_score

# Save the AnnData object to a .h5ad file
adata$write(output_file)

# CTann: Cell Type Annotation using Azimuth

This repository contains the implementation of a **cell type annotation** pipeline using the **Azimuth** tool on single-nucleus RNA sequencing (snRNA-seq) data. The process includes running Azimuth for cell type annotation, analyzing the results, and visualizing the predicted cell types at different levels.

## Dataset

The dataset used for this analysis is **10X snRNA-seq data of the human kidney** (left). It is part of the **HuBMAP** consortium and is available [here](https://portal.hubmapconsortium.org/browse/dataset/180dd06ae35a0fc771f0afe5deefcf23?redirected=True&redirectedFromId=HBM593.CLXN.573&redirectedFromPipeline=Salmon#section-salmon-published).

- **Dataset**: HBM572.JLDL.664
- **Type**: snRNA-seq (10x Genomics v3)
- **Organ**: Kidney (Left)
- **Group**: University of California San Diego TMC
- **Consortium**: HuBMAP

### **Important Note**:
- The data file **`expr.h5ad`** is too large to be uploaded directly to this repository (more than 100 MB).
- To run the analysis, you **must download the dataset manually** from the link above and place it in the `data` folder in this repository.
    - The expected file path is: `data/expr.h5ad`.

## Overview

The main steps involved in this analysis are:

1. **Installation of Required Dependencies**:
    - **Python Libraries**:
        - `anndata`, `requests`, `matplotlib`, `seaborn`, etc. (for data manipulation and visualization).
    - **R Libraries**:
        - `Seurat`, `Azimuth`, `SeuratData`, `reticulate`, etc. (for running the Azimuth tool).

2. **Data Acquisition**:
    - The dataset is downloaded from the [HuBMAP Data Portal](https://portal.hubmapconsortium.org/browse/dataset/180dd06ae35a0fc771f0afe5deefcf23?redirected=True&redirectedFromId=HBM593.CLXN.573&redirectedFromPipeline=Salmon#section-salmon-published).
    - The `expr.h5ad` file, containing the gene expression matrix, is used as input for cell type annotation. **Note**: Since the dataset is larger than 100 MB, it is not uploaded to this repository. You will need to download the dataset manually from the provided link and place it in the `data` folder.

3. **Running Azimuth for Cell Type Annotation**:
    - The Azimuth tool is used to perform cell type annotation at multiple levels (`l1`, `l2`, and `l3`).
    - A custom R script is executed that runs the Azimuth tool and produces predictions for each cell in the dataset.
    - The results are saved into a new `.h5ad` file, which is then read into Python using the `anndata` library.

4. **Data Analysis**:
    - The annotated `.h5ad` file is read into Python, and cell type predictions at various levels (`l1`, `l2`, `l3`) are examined.
    - The distribution of predicted cell types is analyzed, and the number of cells per cell type is computed for each level.

5. **Data Visualization**:
    - **UMAP Plots** are created to visualize the clustering of cells according to predicted cell types at levels `l1`, `l2`, and `l3`.
    - **Hierarchical Clustering** of cells is performed to observe the relationships between different clusters based on their gene expression profiles.
    - **Distribution of Mapping Scores** is visualized to assess the quality of the predictions and the alignment of cells with the reference dataset.

## Running the Notebook

### 1. Install Python and R Dependencies

- **Install Python Libraries**:
    Run the following command to install all required Python libraries:
    ```bash
    pip install anndata requests matplotlib pandas seaborn rpy2
    ```

- **Install R**:
    To run the Azimuth tool, you need to have R installed on your system. You can install R from [here](https://cran.r-project.org/).

- **Install R Libraries**:
    After installing R, install the necessary R libraries. You can use `devtools` to install the required packages directly from GitHub:
    ```R
    # Install Seurat, SeuratData, and Azimuth for Seurat v5
    devtools::install_github("satijalab/seurat", "seurat5")
    devtools::install_github("satijalab/seurat-data", "seurat5")
    devtools::install_github("satijalab/azimuth", "seurat5")
    ```

- **Install R Python Integration**:
    Install the `reticulate` package in R to enable Python-R integration:
    ```R
    install.packages("reticulate")
    ```

### 2. Execute the Jupyter Notebook

- The notebook `azimuth_cell_type_annotation.ipynb` contains the full pipeline. Simply run the notebook from top to bottom in a Jupyter environment.

### 3. R Script Execution

- The notebook automatically generates an R script (`run_azimuth.R`) for running the Azimuth tool. It will execute the script to perform the cell type annotation.

### 4. Data Analysis & Visualization

- The notebook includes detailed steps to analyze the annotated data and visualize the results, including UMAP plots and distribution analysis.

## Outputs

- **Annotated `.h5ad` File**: The final annotated dataset is saved as `annotated_expr.h5ad` which includes predictions for cell types at different levels.
- **UMAP Plots**: UMAP visualizations for predicted cell types at levels `l1`, `l2`, and `l3`.
- **Clustering Analysis**: Hierarchical clustering results for the cells.
- **Mapping Score Distribution**: Histogram displaying the distribution of mapping scores.

## Conclusion

This notebook provides a complete pipeline for performing cell type annotation using the **Azimuth** tool, reading the resulting annotated data, and visualizing the cell types at different levels of granularity. This analysis is useful for exploring the cellular composition of tissues based on single-nucleus RNA sequencing data.

library(Seurat)

# Assuming 'seurat_object' is your Seurat object with cells labeled as 'N', 'T', 'AT', 'NAT'

# Find differentially expressed genes between T cells and other cell types
markers_T_vs_N <- FindMarkers(seurat_object, ident.1 = "T", ident.2 = "N")
markers_T_vs_AT <- FindMarkers(seurat_object, ident.1 = "T", ident.2 = "AT")
markers_T_vs_NAT <- FindMarkers(seurat_object, ident.1 = "T", ident.2 = "NAT")

# Combine the markers into a single list (you may want to adjust this based on your specific criteria)
combined_markers <- list(markers_T_vs_N, markers_T_vs_AT, markers_T_vs_NAT)

# Now you can check the expression of these markers in your AT and NAT cells

#####################

library(Seurat)

# Assuming 'seurat_object' is your Seurat object with cells labeled as 'N', 'T', 'AT', 'NAT'

# Calculate the average expression profile for T cells
t_cells <- subset(seurat_object, idents = "T")
t_cells_avg <- AverageExpression(t_cells)

# Normalize and scale the data if not already done
seurat_object <- NormalizeData(seurat_object)
seurat_object <- ScaleData(seurat_object)

# Extract the scaled data for similarity calculations
data_matrix <- GetAssayData(seurat_object, assay = "RNA", slot = "scale.data")

# Function to calculate cosine similarity
cosine_similarity <- function(x, y) {
  sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
}

# Calculate similarity scores for AT and NAT cells
similarity_scores <- apply(data_matrix, 2, function(cell_expression) {
  cosine_similarity(cell_expression, t_cells_avg$RNA[,1])
})

# Add similarity scores to the metadata of the Seurat object
seurat_object[["cosine_similarity"]] <- similarity_scores

# Now you can access the similarity scores for each cell
similarity_scores_AT <- seurat_object[["cosine_similarity"]][seurat_object$cell_type == "AT"]
similarity_scores_NAT <- seurat_object[["cosine_similarity"]][seurat_object$cell_type == "NAT"]

# You can set a threshold to determine which cells are considered similar to T cells
threshold <- 0.5 # This is an arbitrary example, you'll need to determine an appropriate threshold
similar_cells_AT <- names(similarity_scores_AT[similarity_scores_AT > threshold])
similar_cells_NAT <- names(similarity_scores_NAT[similarity_scores_NAT > threshold])

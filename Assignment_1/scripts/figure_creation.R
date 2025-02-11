library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(reshape2)

# PCA PLOT ----
# Set file paths for eigenvector file and popmap file for creation of PCA plot
eigenvector_file <- "../data/pca_results.eigenvec"
popmap_file <- "../data/popmap.tsv"

eigenvectors <- read.table(eigenvector_file, header = FALSE)

#Adjust the naming conventions for the eigenvectors file
eigenvectors_filtered <- eigenvectors %>% 
  select(!V1) %>% 
  rename(
    Sample = V2, 
    PC1 = V3, 
    PC2 = V4,
    PC3 = V5,
    PC4 = V6,
    PC5 = V7,
    PC6 = V8,
    PC7 = V9,
    PC8 = V10,
    PC9 = V11,
    PC10 = V12
  )

# Read the popmap file and adjust naming conventions
popmap <- read.table(popmap_file, header = FALSE, sep = "\t")
colnames(popmap) <- c("Sample", "Population")

# Merge eigenvectors and popmap by the "Sample" column for creation of PCA plot
pca_data <- merge(eigenvectors_filtered, popmap, by = "Sample")


# Plot PCA 
ggplot(pca_data, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(size = 4) +  # Adjust size of the dots
  labs(title = "Principal Component Analysis of Samples by Population", x = "PC1", y = "PC2") +  # Updated title
  theme_minimal() +
  theme(
    legend.text = element_text(size = 14),   # Adjust size of legend text
    legend.title = element_text(size = 16),  # Adjust size of legend title
    legend.key.size = unit(1.5, "cm"),       # Adjust size of legend keys
    axis.title.x = element_text(size = 16),  # Adjust size of x-axis title
    axis.title.y = element_text(size = 16),  # Adjust size of y-axis title
    plot.title = element_text(size = 20)     # Increase main title size
  )



# PLOT F-STATS DATA IN HEATMAP ----
# Load the F-statistics data
fst_data <- read.table("../data/populations.fst_summary.tsv", header = TRUE, sep = "\t")

# reshape the fst_data for easier retrieval of parameters for matrix
fst_long <- fst_data %>%
  pivot_longer(cols = -X, names_to = "comparison", values_to = "Fst") %>%
  drop_na()

# Create a vector of population names
populations <- c("cs", "pcr", "sj", "wc")

# Create an empty 4x4 matrix (initializing with NA)
fst_matrix_full <- matrix(NA, nrow = 4, ncol = 4)
rownames(fst_matrix_full) <- populations
colnames(fst_matrix_full) <- populations

# Iterate over the data in fst_long to create a matrix for visualization via a heatmap
for(i in 1:nrow(fst_long)) {
  pop1 <- fst_long$X[i]
  pop2 <- fst_long$comparison[i]
  fst_value <- fst_long$Fst[i]
  
  # Fill the matrix with the Fst value for each pair
  fst_matrix_full[pop1, pop2] <- fst_value
  fst_matrix_full[pop2, pop1] <- fst_value
}


# Create the final fst matrix for PCA plot
final_fst_matrix <- as.data.frame(as.table(fst_matrix_full))

# Rename the columns for clarity
colnames(final_fst_matrix) <- c("Population1", "Population2", "Fst")

# Plot the heatmap 
ggplot(final_fst_matrix, aes(x = Population1, y = Population2, fill = Fst)) +
  geom_tile(color = "black") +  
  scale_fill_gradient(low = "#F0E5DE", high = "#C44D58") + 
  labs(title = "Heatmap of Fst Values Between Populations", fill = "Fst") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),  
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 14, face = "bold"),  
        axis.text.y = element_text(size = 14, face = "bold"),  
        legend.title = element_text(size = 14, face = "bold"), 
        legend.text = element_text(size = 12), 
        plot.title = element_text(size = 27, face = "bold", hjust = 0.5)) 



